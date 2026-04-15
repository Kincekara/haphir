version 1.0

task combine_asms {
    input {
        String id
        File flye_asm
        File hifiasm_asm
        File wtdbg2_asm
        File raven_asm
    }

    command <<<
        set -euo pipefail

        # version
        autocycler --version | cut -d " " -f2 > VERSION

        # collect assemblies
        mkdir assemblies
        cp ~{hifiasm_asm} ~{flye_asm} ~{wtdbg2_asm} ~{raven_asm} assemblies/
        # compress
        autocycler compress -i assemblies -a autocycler_out
        # cluster
        autocycler cluster -a autocycler_out
        # trim and resolve
        for c in autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c "$c"
            autocycler resolve -c "$c"
        done
        # combine
        autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

        # rename outputs
        mv autocycler_out/consensus_assembly.fasta ~{id}.autocycler.fasta
        mv autocycler_out/consensus_assembly.gfa ~{id}.autocycler.gfa
    >>>

    output {
        String autocycler_version = read_string("VERSION")
        File assembly_fasta = "~{id}.autocycler.fasta"
        File assembly_graph = "~{id}.autocycler.gfa"
    }

    runtime {
        docker: "staphb/autocycler:0.6.1"
        cpu: 8
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}

