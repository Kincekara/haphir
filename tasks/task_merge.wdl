version 1.0

task merge_asms {
    input {
        String id
        File autocycler_asm
        File plassembler_asm
        File overlaps_paf
        Float identity_threshold = 0.95
        Float coverage_threshold = 0.60
    }

    command <<<
    cat << 'EOF' > merge_plasmids.py
    import argparse
    from collections import defaultdict
    from Bio import SeqIO

    # Thresholds for alignment quality checks
    ID_THRESH = 0.95      # Minimum sequence identity (95%) for considering sequences as matches
    COV_THRESH = 0.60     # Minimum alignment coverage (60%) required for accepting a match

    def parse_header_meta(description):
        """Extract key=value metadata pairs from FASTA header description.
        
        Args:
            description: Full FASTA header description string
            
        Returns:
            Dictionary of metadata key-value pairs found in header
        """
        parts = description.split()
        meta = {}
        # Skip first part (sequence ID), parse remaining parts as key=value pairs
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)  # Split on first '=' only to handle values with '='
                meta[k] = v
        return meta

    def is_circular(meta):
        """Check if a sequence is marked as circular in its metadata.
        
        Args:
            meta: Metadata dictionary (from parse_header_meta)
            
        Returns:
            Boolean indicating if sequence is circular
        """
        return meta.get("circular", "").lower() == "true"

    def read_fasta_with_meta(path):
        """Read FASTA file and extract both sequence records and metadata.
        
        Args:
            path: Path to FASTA file
            
        Returns:
            Tuple of (records dict, metadata dict) indexed by sequence ID
        """
        records = {}  # Store SeqRecord objects
        meta = {}     # Store parsed metadata from headers
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = rec
            meta[rec.id] = parse_header_meta(rec.description)
        return records, meta



    def parse_paf_grouped(paf_path):
        """Parse PAF file and group alignments by (query, target) pairs.
        
        Aggregates multiple alignments between the same query-target pair by:
        - Summing alignment lengths and matching bases
        - Computing combined identity (total_matches / total_alignment_length)
        - Computing combined coverage (total_alignment_length / shorter_sequence_length)
        
        Args:
            paf_path: Path to PAF file containing alignment results
            
        Returns:
            Dictionary mapping (query, target) tuples to aggregated alignment info
        """
        grouped = {}  # Key: (qname, tname), Value: aggregated stats
        
        with open(paf_path) as f:
            for line in f:
                if not line.strip():
                    continue
                # PAF format columns (tab-separated)
                cols = line.rstrip("\n").split("\t")
                qname = cols[0]      # Query sequence name
                qlen = int(cols[1])  # Query sequence length
                tname = cols[5]      # Target sequence name
                tlen = int(cols[6])  # Target sequence length
                matches = int(cols[9])  # Number of matching bases
                aln_len = int(cols[10]) # Alignment block length

                if aln_len == 0:
                    continue

                key = (qname, tname)
                
                if key not in grouped:
                    grouped[key] = {
                        "qname": qname,
                        "tname": tname,
                        "qlen": qlen,
                        "tlen": tlen,
                        "total_aln_len": 0,
                        "total_matches": 0,
                    }
                
                grouped[key]["total_aln_len"] += aln_len
                grouped[key]["total_matches"] += matches
        
        # Compute combined identity and coverage for each group
        results = {}
        for key, stats in grouped.items():
            total_aln = stats["total_aln_len"]
            total_match = stats["total_matches"]
            
            if total_aln > 0:
                combined_ident = total_match / total_aln
                combined_cov = total_aln / min(stats["qlen"], stats["tlen"])
            else:
                combined_ident = 0
                combined_cov = 0
            
            results[key] = {
                "query": stats["qname"],
                "target": stats["tname"],
                "query_len": stats["qlen"],
                "target_len": stats["tlen"],
                "total_alignment_length": total_aln,
                "total_matches": total_match,
                "combined_identity": combined_ident,
                "combined_coverage": combined_cov,
            }
        
        return results




    def choose_best_match(pl_id, hits, auto_meta, plas_meta, id_thresh, cov_thresh):
        """Choose which assembly (plassembler or autocycler) to use for overlapping plasmids.
        
        Decision logic (priority order):
        1. Filter hits by identity and coverage thresholds
        2. If circular vs non-circular: prefer circular
        3. If both circular or both non-circular: prefer longer sequence
        4. If same length: prefer autocycler
        
        Args:
            pl_id: Plassembler plasmid ID
            hits: List of alignment hits for this plasmid
            auto_meta: Autocycler sequence metadata
            plas_meta: Plassembler sequence metadata
            id_thresh: Minimum identity threshold
            cov_thresh: Minimum coverage threshold
            
        Returns:
            Tuple of (decision, best_hit, meets_thresh) where decision is a string indicating which
            assembly to use, best_hit is alignment info (even if thresholds not met), and meets_thresh
            is boolean indicating if best_hit meets quality thresholds
        """
        # Find best hit regardless of thresholds (for reporting)
        if not hits:
            return None, None, False
        
        hits_sorted = sorted(hits, key=lambda h: (h["cov"], h["ident"]), reverse=True)
        best_overall = hits_sorted[0]
        
        # Filter hits that meet quality thresholds
        candidates = [h for h in hits if h["ident"] >= id_thresh and h["cov"] >= cov_thresh]
        if not candidates:
            return None, best_overall, False

        # Sort by coverage first (primary), then identity (secondary), in descending order
        candidates.sort(key=lambda h: (h["cov"], h["ident"]), reverse=True)
        best = candidates[0]  # Take best match
        t = best["target"]    # Autocycler sequence ID

        # Extract metadata for decision making
        plas_m = plas_meta[pl_id]  # Plassembler metadata
        auto_m = auto_meta.get(t, {})  # Autocycler metadata

        # Check circularity status
        plas_circ = is_circular(plas_m)
        auto_circ = is_circular(auto_m)

        # Extract sequence lengths
        plas_len = best["qlen"]   # Plassembler sequence length
        auto_len = best["tlen"]   # Autocycler sequence length

        # Decision 1: If only one is circular, prefer the circular one
        if plas_circ and not auto_circ:
            return ("plassembler", best, True)
        if auto_circ and not plas_circ:
            return ("autocycler", best, True)

        # Decision 2: Both circular or both non-circular - prefer longer sequence
        if plas_len > auto_len:
            return ("plassembler", best, True)
        if auto_len > plas_len:
            return ("autocycler", best, True)

        # Decision 3: Same length - prefer autocycler by default
        return ("autocycler_same_len", best, True)

    def grouped_to_hits(grouped_alignments, query_id):
        """Convert grouped alignment data to hits format for a specific query.
        
        Extracts all (query, target) pairs for a given query ID and converts them
        to the hit format expected by choose_best_match().
        
        Args:
            grouped_alignments: Dictionary from parse_paf_grouped()
            query_id: The query sequence ID to extract hits for
            
        Returns:
            List of hit dictionaries with keys: target, qlen, tlen, ident, cov
        """
        hits = []
        for (query, target), stats in grouped_alignments.items():
            if query == query_id:
                hits.append({
                    "target": target,
                    "qlen": stats["query_len"],
                    "tlen": stats["target_len"],
                    "ident": stats["combined_identity"],
                    "cov": stats["combined_coverage"],
                })
        return hits

    def main():
        ap = argparse.ArgumentParser(
            description="Merge plasmids from plassembler and autocycler, grouping alignments by (query, target) pairs",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
            Example:
                %(prog)s --autocycler auto.fasta --plassembler plas.fasta \\
                        --paf alignments.paf --out merged.fasta
            """
        )
        ap.add_argument("--autocycler", required=True,
                        help="Autocycler FASTA file")
        ap.add_argument("--plassembler", required=True,
                        help="Plassembler FASTA file")
        ap.add_argument("--paf", required=True,
                        help="PAF alignment file")
        ap.add_argument("--out", required=True,
                        help="Output merged FASTA file")
        ap.add_argument("--identity", type=float, default=ID_THRESH,
                        help=f"minimum alignment identity threshold (default: {ID_THRESH})")
        ap.add_argument("--coverage", type=float, default=COV_THRESH,
                        help=f"minimum alignment coverage threshold (default: {COV_THRESH})")
        args = ap.parse_args()

        # Import BioPython for FASTA handling
        from Bio import SeqIO
        
        # Parse grouped PAF alignments
        grouped = parse_paf_grouped(args.paf)

        # Read input files: autocycler and plassembler assemblies
        auto_recs, auto_meta = read_fasta_with_meta(args.autocycler)
        plas_recs, plas_meta = read_fasta_with_meta(args.plassembler)

        # Track which sequences to include in final output
        chosen_autocycler = set(auto_recs.keys())  # Start with all autocycler sequences
        autocycler_choices = {}  # Map autocycler ID to (chosen_plassembler_ID, best_alignment_info)
        extra_plassembler = {}  # Map plassembler ID to (best_match_info) for candidates to include
        updated_auto_meta = {k: dict(v) for k, v in auto_meta.items()}  # Copy of autocycler metadata to potentially update

        # Track decisions for summary report
        decision_rows = []

        # Process each plassembler plasmid to decide inclusion in final output
        for pl_id, pl_rec in plas_recs.items():
            # Find all alignments for this plassembler plasmid against autocycler sequences
            pl_hits = grouped_to_hits(grouped, pl_id)
            decision, best, meets_thresh = choose_best_match(pl_id, pl_hits, auto_meta, plas_meta, args.identity, args.coverage)

            # Handle case where plassembler plasmid doesn't match any autocycler sequence
            if decision is None:
                # best will contain the best hit even if it doesn't meet thresholds
                if best:
                    ident_str = f"{best['ident']:.4f}"
                    cov_str = f"{best['cov']:.4f}"
                    target = best["target"]
                else:
                    ident_str = "-"
                    cov_str = "-"
                    target = "-"
                
                # Keep circular plasmids without matches (likely real plasmids missed by autocycler)
                if is_circular(plas_meta[pl_id]):
                    extra_plassembler[pl_id] = None
                    decision_rows.append([pl_id, target, ident_str, cov_str, "kept_plassembler","no similar hit; circular plasmid"])
                else:
                    # Discard non-circular plasmids without matches (likely spurious)
                    decision_rows.append([pl_id, target, ident_str, cov_str, "ignored","no similar hit; non-circular"])
                continue

            # Get the matching autocycler sequence
            target = best["target"]

            # Case 1: Plassembler version is better (circular and autocycler isn't, or longer)
            if decision == "plassembler":
                # Compare with existing choice for same autocycler, keep the best match
                if target in autocycler_choices:
                    # Get the previous candidate's alignment info
                    prev_pl_id, prev_best = autocycler_choices[target]
                    # Compare by coverage first (primary), then identity (secondary)
                    if (best["cov"], best["ident"]) > (prev_best["cov"], prev_best["ident"]):
                        # New candidate is better, replace previous choice
                        del extra_plassembler[prev_pl_id]
                        autocycler_choices[target] = (pl_id, best)
                        extra_plassembler[pl_id] = best
                    # Otherwise keep previous choice (it's better or equal)
                else:
                    # First plassembler to claim this autocycler
                    autocycler_choices[target] = (pl_id, best)
                    extra_plassembler[pl_id] = best
                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "plassembler", "plassembler preferred (circular or longer)"
                ])

            # Case 2: Autocycler version is better or equal (use it, ignore plassembler)
            elif decision == "autocycler":
                plas_m = plas_meta[pl_id]

                # Ensure target metadata exists before updating (safety check)
                if target not in updated_auto_meta:
                    updated_auto_meta[target] = {}

                # Merge useful metadata from plassembler into autocycler version
                # Copy both plasmid copy-number estimates from plassembler
                for key in ["plasmid_copy_number_short", "plasmid_copy_number_long"]:
                    if key in plas_m:
                        updated_auto_meta[target][key] = plas_m[key]

                # Keep autocycler version, discard plassembler
                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "autocycler", "autocycler preferred (circular or longer); copied short/long copy numbers"
                ])

            # Case 3: Same length sequences - keep autocycler but copy plasmid copy-number metadata
            elif decision == "autocycler_same_len":
                plas_m = plas_meta[pl_id]

                # Ensure target metadata exists before updating (safety check)
                if target not in updated_auto_meta:
                    updated_auto_meta[target] = {}

                # Merge useful metadata from plassembler into autocycler version
                # Copy both plasmid copy-number estimates from plassembler
                for key in ["plasmid_copy_number_short", "plasmid_copy_number_long"]:
                    if key in plas_m:
                        updated_auto_meta[target][key] = plas_m[key]

                decision_rows.append([
                    pl_id, target, f"{best['ident']:.4f}", f"{best['cov']:.4f}",
                    "autocycler_same_len",
                    "same length; autocycler kept; copied short/long copy numbers"
                ])

            # Error case: Unexpected decision value (should never happen)
            else:
                raise ValueError(f"Unexpected decision value for {pl_id}: {decision}")
        
        # Remove autocycler sequences that were chosen as plassembler replacements
        for target, (chosen_pl_id, best_align) in autocycler_choices.items():
            if target in chosen_autocycler:
                chosen_autocycler.remove(target)

        # ===== Build final merged assembly =====
        final_records = []  # Accumulate all sequences for final output

        # Add chosen Autocycler contigs (with possibly updated metadata from plassembler)
        for aid in chosen_autocycler:
            rec = auto_recs[aid]
            meta = updated_auto_meta.get(aid, {})
            meta_str = " ".join(f"{k}={v}" for k, v in meta.items())
            rec.description = f"{rec.id} {meta_str}".strip()
            final_records.append(rec)

        # Add extra Plassembler plasmids not present in Autocycler
        for pid in extra_plassembler:
            rec = plas_recs[pid]
            meta = plas_meta[pid]
            meta_str = " ".join(f"{k}={v}" for k, v in meta.items())
            rec.description = f"{rec.id} {meta_str}".strip()
            final_records.append(rec)

        # Sort all sequences by length (descending) - largest first
        final_records.sort(key=lambda r: len(r.seq), reverse=True)

        # Rename sequences for NCBI submission format
        # First (longest) = chromosome; rest = plasmid1, plasmid2, ...
        chromosome_assigned = False
        plasmid_counter = 1

        for rec in final_records:
            # Extract and preserve metadata from description
            parts = rec.description.split()
            meta = " ".join(parts[1:])  # Keep everything after original ID

            if not chromosome_assigned:
                # Longest sequence is the chromosome
                rec.id = "chromosome"
                rec.description = f"chromosome {meta}".strip()
                chromosome_assigned = True
            else:
                # Remaining sequences are numbered plasmids
                rec.id = f"plasmid{plasmid_counter}"
                rec.description = f"{rec.id} {meta}".strip()
                plasmid_counter += 1

        # Write final merged FASTA to output file
        SeqIO.write(final_records, args.out, "fasta")

        # Write decision summary as TSV for debugging and analysis
        with open("merge_summary.tsv", "w") as out:
            # Header: plasmid ID, matching sequence, alignment metrics, decision, reasoning
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
    --out ~{id}.merged.fasta \
    --identity ~{identity_threshold} \
    --coverage ~{coverage_threshold}

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

