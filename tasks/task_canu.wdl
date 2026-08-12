version 1.0

task hicanu_asm {
    input {
        String id
        String genome_size
        File long_fq
        Int cpu = 8
    }

    command <<<
        set -euo pipefail

        # version 
        canu -version | cut -d ' ' -f 2 > VERSION
    
        # assemble with hicanu
        canu \
        -p ~{id} \
        -d out \
        genomeSize=~{genome_size} \
        -pacbio-hifi ~{long_fq} \
        useGrid=false \
        maxThreads=~{cpu}

        # get contig lengths
        mv out/~{id}.contigs.fasta ~{id}.hicanu.fasta
        echo "HiCanu" > ~{id}.hicanu.ctg_len.txt
        awk -F"len=| " '/^>/{print $3}' ~{id}.hicanu.fasta | sort -nr >> ~{id}.hicanu.ctg_len.txt
    >>>

    output {
        String canu_version = read_string("VERSION")
        File assembly_fasta = "~{id}.hicanu.fasta"
        File ctg_len = "~{id}.hicanu.ctg_len.txt"
    }

    runtime {
        docker: "staphb/canu:2.3"
        cpu: cpu
        memory: "32 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}