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
        # give contigs from Hifiasm and Flye extra consensus weight
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' assemblies/*hifiasm.fasta
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' assemblies/*flye.fasta
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
        
        # get contig lengths
        echo "Autocycler" > ~{id}.autocycler.ctg_len.txt
        awk -F'length=' '/^>/{split($2,a," "); print a[1]}' ~{id}.autocycler.fasta | sort -nr >> ~{id}.autocycler.ctg_len.txt
    >>>

    output {
        String autocycler_version = read_string("VERSION")
        File assembly_fasta = "~{id}.autocycler.fasta"
        File assembly_graph = "~{id}.autocycler.gfa"
        File ctg_len = "~{id}.autocycler.ctg_len.txt"
    }

    runtime {
        docker: "staphb/autocycler:0.6.2"
        cpu: 8
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}

