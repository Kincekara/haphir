version 1.0

task merge_asms {
    input {
        String id
        File autocycler_asm
        File plassembler_asm
        File overlaps_paf
    }

    command <<<
    cat << 'EOF' > merge_plasmids.py
    import argparse
    from collections import defaultdict
    from Bio import SeqIO

    ID_THRESH = 0.95
    COV_THRESH = 0.60

    def parse_header_meta(description):
        parts = description.split()
        meta = {}
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                meta[k] = v
        return meta

    def is_circular(meta):
        return meta.get("circular", "").lower() == "true"

    def read_fasta_with_meta(path):
        records = {}
        meta = {}
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = rec
            meta[rec.id] = parse_header_meta(rec.description)
        return records, meta

    def parse_paf(paf_path):
        hits = defaultdict(list)
        with open(paf_path) as f:
            for line in f:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                qname = cols[0]
                qlen = int(cols[1])
                tname = cols[5]
                tlen = int(cols[6])
                aln_len = int(cols[9])
                matches = int(cols[10])

                if aln_len == 0:
                    continue

                ident = matches / aln_len
                cov = aln_len / min(qlen, tlen)

                hits[qname].append({
                    "target": tname,
                    "qlen": qlen,
                    "tlen": tlen,
                    "aln_len": aln_len,
                    "ident": ident,
                    "cov": cov,
                })
        return hits

    def choose_best_match(pl_id, hits, auto_meta, plas_meta):
        candidates = [h for h in hits if h["ident"] >= ID_THRESH and h["cov"] >= COV_THRESH]
        if not candidates:
            return None, None

        candidates.sort(key=lambda h: (h["cov"], h["ident"]), reverse=True)
        best = candidates[0]
        t = best["target"]

        plas_m = plas_meta[pl_id]
        auto_m = auto_meta.get(t, {})

        plas_circ = is_circular(plas_m)
        auto_circ = is_circular(auto_m)

        plas_len = best["qlen"]
        auto_len = best["tlen"]

        if plas_circ and not auto_circ:
            return ("plassembler", best)
        if auto_circ and not plas_circ:
            return ("autocycler", best)

        if plas_len > auto_len:
            return ("plassembler", best)
        if auto_len > plas_len:
            return ("autocycler", best)

        return ("autocycler_same_len", best)

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--autocycler", required=True)
        ap.add_argument("--plassembler", required=True)
        ap.add_argument("--paf", required=True)
        ap.add_argument("--out", required=True)
        args = ap.parse_args()

        auto_recs, auto_meta = read_fasta_with_meta(args.autocycler)
        plas_recs, plas_meta = read_fasta_with_meta(args.plassembler)
        hits = parse_paf(args.paf)

        chosen_autocycler = set(auto_recs.keys())
        extra_plassembler = []
        updated_auto_meta = {k: dict(v) for k, v in auto_meta.items()}

        decision_rows = []

        for pl_id, pl_rec in plas_recs.items():
            pl_hits = hits.get(pl_id, [])
            decision, best = choose_best_match(pl_id, pl_hits, auto_meta, plas_meta)

            if decision is None:
                if is_circular(plas_meta[pl_id]):
                    extra_plassembler.append(pl_id)
                    decision_rows.append([pl_id, "-", "-", "-", "kept_plassembler","no similar hit; circular plasmid"])
                else:
                    decision_rows.append([pl_id, "-", "-", "-", "ignored","no similar hit; non-circular"])
                continue

            target = best["target"]

            if decision == "plassembler":
                if target in chosen_autocycler:
                    chosen_autocycler.remove(target)
                extra_plassembler.append(pl_id)
                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "plassembler", "plassembler preferred (circular or longer)"
                ])
                continue

            if decision == "autocycler":
                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "autocycler", "autocycler preferred (circular or longer)"
                ])
                continue

            elif decision == "autocycler_same_len":
                plas_m = plas_meta[pl_id]

                # Copy both Plassembler copy-number fields if present
                for key in ["plasmid_copy_number_short", "plasmid_copy_number_long"]:
                    if key in plas_m:
                        updated_auto_meta[target][key] = plas_m[key]

                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "autocycler_same_len",
                    "same length; autocycler kept; copied short/long copy numbers"
                ])
                continue

        # Build final records
        final_records = []

        # Autocycler contigs (with possibly updated metadata)
        for aid in chosen_autocycler:
            rec = auto_recs[aid]
            meta = updated_auto_meta.get(aid, {})
            meta_str = " ".join(f"{k}={v}" for k, v in meta.items())
            rec.description = f"{rec.id} {meta_str}".strip()
            final_records.append(rec)

        # Extra Plassembler plasmids
        for pid in extra_plassembler:
            rec = plas_recs[pid]
            meta = plas_meta[pid]
            meta_str = " ".join(f"{k}={v}" for k, v in meta.items())
            rec.description = f"{rec.id} {meta_str}".strip()
            final_records.append(rec)

        # Sort contigs by length (descending)
        final_records.sort(key=lambda r: len(r.seq), reverse=True)

        # Rename headers for NCBI submission
        chromosome_assigned = False
        plasmid_counter = 1

        for rec in final_records:
            # Extract metadata from description
            parts = rec.description.split()
            meta = " ".join(parts[1:])  # everything after original ID

            if not chromosome_assigned:
                # First (longest) contig → chromosome
                rec.id = "chromosome"
                rec.description = f"chromosome {meta}".strip()
                chromosome_assigned = True
            else:
                # Remaining contigs → plasmid1, plasmid2, ...
                rec.id = f"plasmid{plasmid_counter}"
                rec.description = f"{rec.id} {meta}".strip()
                plasmid_counter += 1

        # Write final FASTA
        SeqIO.write(final_records, args.out, "fasta")

        # Write TSV summary
        with open("merge_summary.tsv", "w") as out:
            out.write("plasmid_id\tbest_hit\tidentity\tcoverage\tdecision\treason\n")
            for row in decision_rows:
                out.write("\t".join(row) + "\n")

    if __name__ == "__main__":
        main()
    EOF

        chmod +x merge_plasmids.py

        python3 merge_plasmids.py \
        --autocycler ~{autocycler_asm} \
        --plassembler ~{plassembler_asm} \
        --paf ~{overlaps_paf} \
        --out ~{id}.merged.fasta

        mv merge_summary.tsv ~{id}.merge_summary.tsv
    >>>

    output {
        File merged_fasta = "~{id}.merged.fasta"
        File merge_summary = "~{id}.merge_summary.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/biopython:1.84"
        cpu: 1
        memory: "4 GiB"
    }


}

