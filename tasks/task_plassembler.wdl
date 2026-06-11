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
        if [ -s "~{short_fq1}" ] && [ -s "~{short_fq2}" ]; then
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
        else
            plassembler long \
            --threads ~{cpu} \
            --database /plassembler_db \
            --pacbio_model pacbio-hifi \
            --longreads ~{long_fq} \
            --flye_assembly ~{flye_asm} \
            --flye_info ~{flye_info} \
            --skip_qc \
            --prefix ~{id} \
            --outdir out
        fi

        # get contig lengths
        echo "Plassembler" > ~{id}.plassembler.ctg_len.txt
        awk 'NR > 1 {print $2}' out/~{id}_summary.tsv >> ~{id}.plassembler.ctg_len.txt
    >>>

    output {
        String plassembler_version = read_string("VERSION")
        File plasmids = "out/~{id}_plasmids.fasta"
        File graph = "out/~{id}_plasmids.gfa"
        File summary = "out/~{id}_summary.tsv"
        File ctg_len = "~{id}.plassembler.ctg_len.txt"
    }

    runtime {
        docker: "staphb/plassembler:1.8.1"
        cpu: cpu
        memory: "16 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}