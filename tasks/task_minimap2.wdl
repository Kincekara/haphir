version 1.0

task label_and_align {
    input {
        String id
        File autocycler_asm
        File plassembler_asm
    }

    command <<<
        set -euo pipefail

        # function to mark sequence headers with a prefix
        mark_headers() {
            local input_file=$1
            local output_file=$2
            local prefix=$3
            
            awk -v prefix="$prefix" '
            /^>/ {
                count++
                # Remove leading ">" and split header into fields
                sub(/^>/, "", $0)
                split($0, a, " ")
                # Replace only the first field with prefixN
                a[1] = prefix count
                # Reconstruct header
                printf ">%s", a[1]
                for (i=2; i<=length(a); i++) printf " %s", a[i]
                printf "\n"
                next
            }
            {print}
            ' "$input_file" > "$output_file"
        }

        # modify headers
        mark_headers ~{autocycler_asm} ~{id}.autocycler.marked.fasta "autocycler"
        mark_headers ~{plassembler_asm} ~{id}.plasmids.marked.fasta "plassembler"

        # version 
        minimap2 --version > VERSION

        # align with minimap2
        minimap2 -x asm5 ~{id}.autocycler.marked.fasta ~{id}.plasmids.marked.fasta > ~{id}.overlaps.paf
    >>>

    output {
        String minimap_version = read_string("VERSION")
        File overlaps_paf = "~{id}.overlaps.paf"
        File autocycler_fasta = "~{id}.autocycler.marked.fasta"
        File plasmids_fasta = "~{id}.plasmids.marked.fasta"
    }

    runtime {
        docker: "staphb/minimap2:2.30"
        cpu: 2
        memory: "4 GiB"
    }
}
