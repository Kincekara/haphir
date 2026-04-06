version 1.0

task trim_pe {
    input {
        String id
        File? short_fq1
        File? short_fq2
        Int cpu = 4
    }

    command <<<
        set -euo pipefail

        # version 
        fastp --version > VERSION

        # trim reads with fastp
        fastp \
        -i ~{short_fq1} \
        -I ~{short_fq2} \
        -o ~{id}.trimmed.fq1.gz \
        -O ~{id}.trimmed.fq2.gz \
        --length_required 70 \
        --average_qual 30 \
        --cut_front_window_size 1 \
        --cut_front_mean_quality 10 \
        -3 \
        --cut_tail_window_size 1 \
        --cut_tail_mean_quality 10 \
        -r \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --detect_adapter_for_pe \
        --thread ~{cpu} \
        -h ~{id}.fastp.html
    >>>

    output {
        String fastp_version = read_string("VERSION")
        File short_fq1_trimmed  = "~{id}.trimmed.fq1.gz"
        File short_fq2_trimmed  = "~{id}.trimmed.fq2.gz"
        File html_report = "~{id}.fastp.html"
    }

    runtime {
        container: "staphb/fastp:1.1.0"
        cpu: cpu
        memory: "4 GiB"
        preemptible: 0
        maxRetries: 3
    }
}