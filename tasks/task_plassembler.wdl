version 1.0

task plassembler_asm {
    input {
        String id
        File long_fq
        File? short_fq1
        File? short_fq2       
        File flye_asm
        File flye_info
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        plassembler --version | cut -d " " -f3 | tr -d "\n" > VERSION

        # plassembler
        plassembler run \
        --threads ~{cpu} \
        --database /plassembler_db \
        --pacbio_model pacbio-hifi \
        --longreads ~{long_fq} \
        --short_one ~{short_fq1} \
        --short_two ~{short_fq2} \
        --flye_assembly ~{flye_asm} \
        --flye_info ~{flye_info} \
        --skip_qc \
        --prefix ~{id} \
        --outdir out
    >>>

    output {
        String plassembler_version = read_string("VERSION")
        File plasmids = "out/~{id}_plasmids.fasta"
        File graph = "out/~{id}_plasmids.gfa"
        File summary = "out/~{id}_summary.tsv"
    }

    runtime {
        docker: "staphb/plassembler:1.8.2"
        cpu: cpu
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}