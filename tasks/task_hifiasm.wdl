version 1.1

task hifiasm_asm {
    input {
        String id
        String genome_size
        File long_fq
        Int cpu = 16
    }

    command <<<
        set -euo pipefail

        # version 
        hifiasm --version > VERSION

        # round genome size for hifiasm
        gsize=$(printf "%dm\n" $(( ("~{genome_size}"+500000)/1000000 )))
    
        # assemble with hifiasm
        hifiasm \
        -o ~{id} \
        -t ~{cpu} \
        --hg-size "$gsize" \
        ~{long_fq} 2> hifiasm.log

        # gfa to fasta
        mv ~{id}.bp.p_ctg.gfa ~{id}.hifiasm.gfa
        awk '/^S/{print ">"$2;print $3}' ~{id}.hifiasm.gfa > ~{id}.hifiasm.fasta
    >>>

    output {
        String hifiasm_version = read_string("VERSION")
        File assembly_graph = "~{id}.hifiasm.gfa"
        File assembly_fasta = "~{id}.hifiasm.fasta"
    }

    runtime {
        container: "staphb/hifiasm:0.25.0"
        cpu: cpu
        memory: "32 GiB"
        disks: "local-disk 200 SSD"
        preemptible: 2
        maxRetries: 5
    }
}