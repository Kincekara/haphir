version 1.0

task reorient {
    input {
        String id
        File long_asm
        Int cpu = 2
    }

    command <<<
        set -euo pipefail

        # version 
        dnaapler --version | cut -d " " -f3 > VERSION

        # reorient assembly with dnaapler
        dnaapler all \
        -i ~{long_asm} \
        -o out \
        -t ~{cpu}
        # rename output
        mv out/dnaapler_reoriented.fasta ~{id}.dnaapler.fasta
        mv out/dnaapler_all_reorientation_summary.tsv ~{id}.dnaapler_summary.tsv 
    >>>

    output {
        String dnaapler_version = read_string("VERSION")
        File reoriented_fasta = "~{id}.dnaapler.fasta"
        File dnaapler_summary = "~{id}.dnaapler_summary.tsv"
    }

    runtime {
        docker: "staphb/dnaapler:1.3.0"
        cpu: cpu
        memory: "8 GiB"
        disks: "local-disk 50 SSD"
        preemptible: 0
        maxRetries: 3
    }
}