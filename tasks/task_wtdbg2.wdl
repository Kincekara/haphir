version 1.0

task wtdbg2_asm {
    input {
        String id
        File long_fq
        String genome_size
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        wtdbg2 --version | cut -d " " -f2 > VERSION

        # assemble with wtdb2
        wtdbg2 \
        -t ~{cpu} \
        -i ~{long_fq} \
        -o ~{id} \
        -g ~{genome_size} \
        -x ccs \
        -S 2     

        # derive consensus
        wtpoa-cns \
        -t ~{cpu} \
        -i ~{id}.ctg.lay.gz -fo ~{id}.wtdbg2.fasta

        # get contig lengths
        echo "Wtdbg2" > ~{id}.wtdbg2.ctg_len.txt
        awk -F'len=' '/^>/{print $2}' ~{id}.wtdbg2.fasta | sort -nr >> ~{id}.wtdbg2.ctg_len.txt
    >>>

    output {
        String wtdbg2_version = read_string("VERSION")
        File assembly_fasta = "~{id}.wtdbg2.fasta"
        File ctg_len = "~{id}.wtdbg2.ctg_len.txt"
    }

    runtime {
        docker: "staphb/wtdbg2:2.5"
        cpu: cpu
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}

