version 1.0

task polish {
    input {
        String id
        File draft_asm
        File? short_fq1
        File? short_fq2
        Int cpu = 4
    }

    command <<<
        set -euo pipefail

        # version 
        polypolish --version > VERSION

        # index    
        bwa index ~{draft_asm}
        # map
        bwa mem -t ~{cpu} -a ~{draft_asm} ~{short_fq1} > alignments_1.sam
        bwa mem -t ~{cpu} -a ~{draft_asm} ~{short_fq2} > alignments_2.sam
        # filter
        polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
        # polish
        polypolish polish ~{draft_asm} filtered_1.sam filtered_2.sam > ~{id}.polished.fasta
    >>>

    output {
        String version = read_string("VERSION")
        File polished_fasta = "~{id}.polished.fasta"
    }

    runtime {
        docker: "staphb/polypolish:0.6.1-bwa"
        cpu: cpu
        memory: "8 GiB"
        preemptible: 0
        maxRetries: 3
    }
}