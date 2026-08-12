version 1.0

task lja_asm {
    input {
        String id
        File long_fq
        Int cpu = 8
    }

    command <<<
        set -euo pipefail
    
        # assemble with lja
        lja \
        -t ~{cpu} \
        -o out \
        --reads ~{long_fq} > lja.out.txt 2> lja.err.txt

        # rename output
        mv out/assembly.fasta ~{id}.lja.fasta
        mv out/mdbg.gfa ~{id}.lja.gfa

        # get contig lengths
        echo "LJA" > ~{id}.lja.ctg_len.txt
        awk '/^>/ {if (len) print len; len=0; next} {len += length($0)} END {if (len) print len}' ~{id}.lja.fasta | sort -nr >> ~{id}.lja.ctg_len.txt

    >>>

    output {
        String lja_version = "0.2"
        File assembly_fasta = "~{id}.lja.fasta"
        File assembly_graph = "~{id}.lja.gfa"
        File ctg_len = "~{id}.lja.ctg_len.txt"
    }

    runtime {
        docker: "staphb/lja:0.2-bugfix"
        cpu: cpu
        memory: "32 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}