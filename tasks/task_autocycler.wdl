version 1.0

task combine_asms {
    input {
        String id
        File flye_asm
        File hifiasm_asm
        File raven_asm
        File opt_asm
    }

    command <<<
        set -euo pipefail

        # version
        autocycler --version | cut -d " " -f2 > VERSION

        # create helper functions
        run_autocycler() {
            local input_dir="$1" output_dir="$2" log_file="$3"
            autocycler compress -i "$input_dir" -a "$output_dir"
            autocycler cluster -a "$output_dir"
            for c in "$output_dir"/clustering/qc_pass/cluster_*; do
                [ -d "$c" ] || continue
                autocycler trim -c "$c"
                autocycler resolve -c "$c"
            done
            autocycler combine -a "$output_dir" -i "$output_dir"/clustering/qc_pass/cluster_*/5_final.gfa 2>&1 | tee "$log_file"
        }

        finalize_output() {
            local src_prefix="$1" dst_prefix="$2"
            if [ -f "$src_prefix".fasta ]; then
                mv "$src_prefix".fasta "$dst_prefix".fasta
            fi
            if [ -f "$src_prefix".gfa ]; then
                mv "$src_prefix".gfa "$dst_prefix".gfa
            fi
        }

        clean_short_contigs() {
            local input="$1" output="$2"
            awk 'BEGIN {RS=">"; FS="\n"} NR>1 {seq=""; for (i=2; i<=NF; i++) seq=seq$i; if (length(seq) >= 2000) printf ">%s", $0}' "$input" > "$output"
        }

        # collect assemblies
        mkdir assemblies
        cp ~{hifiasm_asm} ~{flye_asm} ~{raven_asm} ~{opt_asm} assemblies/

        # give extra consensus weight to contigs from Hifiasm and Flye 
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' assemblies/*hifiasm.fasta
        sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' assemblies/*flye.fasta
        
        # run autocycler
        run_autocycler "assemblies" "autocycler_out" "combine.log"
        
        # check if consensus assembly is fully resolved
        if grep -q "Consensus assembly is fully resolved" combine.log; then
            echo "Consensus assembly is fully resolved"
            finalize_output "autocycler_out/consensus_assembly" "~{id}.autocycler"
            echo "SUCCESS" > RESULT
        elif grep -q "One or more clusters failed to fully resolve" combine.log; then
            echo "Checking chromosome..."
            # check if chromosome is fully resolved
            if head -n 1 autocycler_out/consensus_assembly.fasta | grep -q "circular"; then
                echo "Chromosome is fully resolved, trying to clean up the assembly..."
                clean_short_contigs "autocycler_out/consensus_assembly.fasta" "autocycler.clean.fasta"              
                # rename outputs
                mv autocycler.clean.fasta ~{id}.autocycler.fasta
                mv autocycler_out/consensus_assembly.gfa ~{id}.autocycler.gfa
                echo "SUCCESS" > RESULT 
            else
                echo "Chromosome is not fully resolved, assembly failed"
                finalize_output "autocycler_out/consensus_assembly" "~{id}.autocycler"
                echo "FAIL" > RESULT       
            fi
        else
            echo "Autocycler failed to produce a consensus assembly"
            echo "FAIL" > RESULT
        fi              

        # get contig lengths
        echo "Autocycler" > ~{id}.autocycler.ctg_len.txt
        awk -F'length=' '/^>/{split($2,a," "); print a[1]}' ~{id}.autocycler.fasta | sort -nr >> ~{id}.autocycler.ctg_len.txt
    >>>

    output {
        String autocycler_version = read_string("VERSION")
        File assembly_fasta = "~{id}.autocycler.fasta"
        File assembly_graph = "~{id}.autocycler.gfa"
        File ctg_len = "~{id}.autocycler.ctg_len.txt"
        String autocycler_result = read_string("RESULT")
    }

    runtime {
        docker: "staphb/autocycler:0.6.2"
        cpu: 8
        memory: "16 GiB"
        preemptible: 2
        maxRetries: 5
    }
}

