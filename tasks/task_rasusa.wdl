version 1.0

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
        docker: "staphb/rasusa:4.0.0"
        cpu: 2
        memory: "2 GiB"
        preemptible: 0
        maxRetries: 3
    }
}

task downsample_pe {
    input {
        String id
        File? short_fq1
        File? short_fq2
        String genome_size
    }

    command <<<
        set -euo pipefail

        # downsample pe reads
        rasusa reads \
        --seed 42 \
        --coverage 110 \
        --genome-size ~{genome_size} \
        -o ~{id}.downsampled.r1.fastq.gz \
        -o ~{id}.downsampled.r2.fastq.gz \
        ~{short_fq1} ~{short_fq2}
    >>>

    output {
        File ds_short_fq1 = "~{id}.downsampled.r1.fastq.gz"
        File ds_short_fq2 = "~{id}.downsampled.r2.fastq.gz"
    }

    runtime {
        docker: "staphb/rasusa:4.0.0"
        cpu: 2
        memory: "2 GiB"
        preemptible: 0
        maxRetries: 3
    }
    
}