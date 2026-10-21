version 1.0

task filter {
    input {
        String id
        File long_fq
        Int min_lenght = 0 # set 0 to disable length filtering, default 1000
        Int bias = 100 # set 100 to disable base content bias correction, default 1
        Int cpu = 4
    }

    command <<<
        # version
        tgsfilter -h | tail -n2 |  grep -oE '[0-9]+\.[0-9]+' > VERSION

        # adapter filter & qc
        tgsfilter \
        -t ~{cpu} \
        -x hifi \
        -l ~{min_lenght} \
        -b ~{bias} \
        -i ~{long_fq} \
        -o ~{id}.filtered.fastq.gz 2>&1 | tee ~{id}.tgsfilter.log

        # rename html report
        mv ~{id}.filtered.html ~{id}.tgsfilter.html
    >>>

    output {
        String tgsfilter_version = read_string("VERSION")
        File filtered_fq = "~{id}.filtered.fastq.gz"
        File tgsfilter_log = "~{id}.tgsfilter.log"
        File html_report = "~{id}.tgsfilter.html"
    }

    runtime {
        docker: "staphb/tgsfilter:1.11"
        cpu: cpu
        memory: "4 GiB"
        preemptible: 0
        maxRetries: 3
    }
}
