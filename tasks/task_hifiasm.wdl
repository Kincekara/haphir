version 1.0

task hifiasm_asm {
    input {
        String id
        String genome_size
        File long_fq
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        hifiasm --version > VERSION
    
        # assemble with hifiasm
        hifiasm \
        -o ~{id} \
        -t ~{cpu} \
        --hg-size ~{genome_size} \
        ~{long_fq} 2> hifiasm.log

        # gfa to fasta
        mv ~{id}.bp.p_ctg.gfa ~{id}.hifiasm.gfa
        awk '/^S/{print ">"$2" length="length($3);print $3}' ~{id}.hifiasm.gfa > ~{id}.hifiasm.fasta

        # get contig lengths
        echo "Hifiasm" > ~{id}.hifiasm.ctg_len.txt
        awk -F'length=' '/^>/{print $2}' ~{id}.hifiasm.fasta | sort -nr >> ~{id}.hifiasm.ctg_len.txt
    >>>

    output {
        String hifiasm_version = read_string("VERSION")
        File assembly_graph = "~{id}.hifiasm.gfa"
        File assembly_fasta = "~{id}.hifiasm.fasta"
        File ctg_len = "~{id}.hifiasm.ctg_len.txt"
    }

    runtime {
        docker: "staphb/hifiasm:0.25.0"
        cpu: cpu
        memory: "32 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}