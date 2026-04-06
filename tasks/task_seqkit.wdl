version 1.1

task merge_asms {
    input {
        String id
        File? long_asm
        File? plassembler_asm
    }

    command <<<
        set -euo pipefail

        # version
        seqkit version > VERSION

        # extract chromosome
        seqkit head -n 1 ~{long_asm} > chromosome.fasta
        awk '
        /^>/ {
            # Remove leading ">"
            sub(/^>/, "", $0)
            # Split header into fields
            split($0, a, " ")
            # Replace only the first field with "chromosome"
            a[1] = "chromosome"
            # Reconstruct header
            printf ">%s", a[1]
            for (i=2; i<=NF; i++) printf " %s", a[i]
            printf "\n"
            next
        }
        {print}
        ' chromosome.fasta > chromosome_renamed.fasta
        
        # extract plasmids
        seqkit grep -n -r -p "circular=true" ~{plassembler_asm} -o plasmids.fasta
        awk '
        /^>/ {
            count++
            # Remove leading ">" and split header into fields
            sub(/^>/, "", $0)
            split($0, a, " ")
            # Replace only the first field with plasmidN
            a[1] = "plasmid" count
            # Reconstruct header
            printf ">%s", a[1]
            for (i=2; i<=length(a); i++) printf " %s", a[i]
            printf "\n"
            next
        }
        {print}
        ' plasmids.fasta > plasmids_renamed.fasta

        # merge
        cat chromosome_renamed.fasta plasmids_renamed.fasta > ~{id}.final_assembly.fasta
        # stats
        seqkit stats ~{id}.final_assembly.fasta > ~{id}.seqkit.stats.txt
    >>>

    output {
        String seqkit_version = read_string("VERSION")
        File final_assembly = "~{id}.final_assembly.fasta"
        File seqkit_stats = "~{id}.seqkit.stats.txt"
    }

    runtime {
        container: "staphb/seqkit:2.13.0"
        cpu: 1
        memory: "1 GiB"
        preemptible: 0
        maxRetries: 3
    }
}