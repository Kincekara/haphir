version 1.0

task bam_to_fastq {
    input {
        String id
        File bam
    }

    command <<<
        set -euo pipefail

        # bam2fastq version
        bam2fastq --version | head -n 1 | cut -d " " -f2 > VERSION

        # index file
        pbindex ~{bam}
        
        # bam to fastq
        bam2fastq -o ~{id}.hifi ~{bam}
    >>>

    output {
        String bam2fastq_version = read_string("VERSION")
        File long_fq = "~{id}.hifi.fastq.gz"
    }

    runtime {
        docker: "staphb/pbtk:3.5.0"
        memory: "4 GiB"
        cpu: 4
        disks: "local-disk 200 SSD"
        preemptible: 0
        maxRetries: 3
    }
}