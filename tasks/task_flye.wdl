version 1.1

task flye_asm {
    input {
        String id
        File long_fq
        String genome_size
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        flye --version > VERSION

        # assemble with flye
        flye \
        --threads ~{cpu} \
        --pacbio-hifi ~{long_fq} \
        --genome-size ~{genome_size} \
        --out-dir out

        # rename outputs
        mv ./out/assembly.fasta ~{id}.flye.fasta
        mv ./out/assembly_graph.gfa ~{id}.flye.gfa
        mv ./out/assembly_info.txt ~{id}.flye_info.txt
    >>>

    output {
        String flye_version = read_string("VERSION")        
        File assembly_fasta = "~{id}.flye.fasta"
        File assembly_graph = "~{id}.flye.gfa"
        File assembly_info = "~{id}.flye_info.txt"
    }

    runtime {
        container: "staphb/flye:2.9.6"
        cpu: cpu
        memory: "16 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}