version 1.0

task annotation {
    input {
        String id
        File final_asm
        String? organism
        Int cpu = 4       
    }

    command <<<
        set -euo pipefail
        
        # version 
        bakta --version | cut -d " " -f2 > VERSION

        # annotate with bakta
        if [ -n "~{organism}" ]; then
            genus=$(echo ~{organism} | cut -d ' ' -f1)
            species=$(echo ~{organism} | cut -d ' ' -f2)

            bakta \
            --threads ~{cpu} \
            --prefix ~{id} \
            --complete \
            --compliant \
            --output bakta \
            --genus "$genus" \
            --species "$species" \
            ~{final_asm}
        else
            bakta \
            --threads ~{cpu}\
            --prefix ~{id} \
            --complete \
            --compliant \
            --output bakta \
            ~{final_asm}
        fi

        # compress outputs
        tar -czvf ~{id}.bakta.tar.gz bakta/
    >>>

    output {
        String bakta_version = read_string("VERSION")
        File bakta_outputs = "~{id}.bakta.tar.gz"
        File bakta_gff = "bakta/~{id}.gff3"
        File bakta_faa = "bakta/~{id}.faa"
        File bakta_fna = "bakta/~{id}.fna"
    }

    runtime {
        docker: "staphb/bakta:1.11.4-6.0-light"
        cpu: cpu
        memory: "8 GiB"
    }
}