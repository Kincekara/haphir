version 1.1

task downsample {
    input {
        String id
        File long_fq
        String genome_size
    }

    command <<<
        set -euo pipefail

        # version 
        rasusa --version > VERSION
        
        # downsample reads
        rasusa reads \
        --seed 42 \
        --coverage 110 \
        --genome-size ~{genome_size} \
        --output ~{id}.downsampled.fastq.gz \
        ~{long_fq}
    >>>

    output {
        String rasusa_version = read_string("VERSION")
        File downsampled_fq = "~{id}.downsampled.fastq.gz"
    }

    runtime {
        container: "staphb/rasusa:3.0.0"
        cpu: 2
        memory: "2 GiB"
        preemptible: 0
        maxRetries: 3
    }
}