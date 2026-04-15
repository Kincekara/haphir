version 1.0

task amr {
    input {
        String id
        String? organism
        File? assembly
        File? bakta_faa
        File? bakta_gff
        Int cpu = 4
    }

    command <<<
        set -euo pipefail

        # version
        amrfinder --version > VERSION
        amrfinder --database_version | grep "Database version" | cut -d " " -f3 > DB_VERSION

        ## curated organisms ##
        # A. baumannii-calcoaceticus species complex
        declare -a abcc=(
            "Acinetobacter baumannii"
            "Acinetobacter calcoaceticus"
            "Acinetobacter lactucae"
            "Acinetobacter nosocomialis"
            "Acinetobacter pittii"
            "Acinetobacter seifertii"
        )
        # Burkholderia cepacia species complex
        declare -a bcc=(
            "Burkholderia aenigmatica"
            "Burkholderia ambifaria"   
            "Burkholderia anthina"   
            "Burkholderia arboris"   
            "Burkholderia catarinensis"   
            "Burkholderia cenocepacia"   
            "Burkholderia cepacia" 
            "Burkholderia cf. cepacia"  
            "Burkholderia contaminans"   
            "Burkholderia diffusa"   
            "Burkholderia dolosa"   
            "Burkholderia lata"   
            "Burkholderia latens"
            "Burkholderia metallica"  
            "Burkholderia multivorans"   
            "Burkholderia orbicola"   
            "Burkholderia paludis"   
            "Burkholderia pseudomultivorans"   
            "Burkholderia puraquae"   
            "Burkholderia pyrrocinia"   
            "Burkholderia semiarida"   
            "Burkholderia seminalis"   
            "Burkholderia sola"   
            "Burkholderia stabilis"   
            "Burkholderia stagnalis"   
            "Burkholderia territorii"   
            "Burkholderia ubonensis"   
            "Burkholderia vietnamiensis" 
        )
        # Burkholderia pseudomallei species complex
        declare -a bpc=(
            "Burkholderia humptydooensis"   
            "Burkholderia mallei"   
            "Burkholderia mayonis"   
            "Burkholderia oklahomensis"   
            "Burkholderia pseudomallei"   
            "Burkholderia savannae"   
            "Burkholderia singularis"   
            "Burkholderia thailandensis"   
        )
        # other species
        declare -a taxa=(   
            "Citrobacter freundii"
            "Clostridioides difficile"
            "Enterobacter asburiae"
            "Enterobacter cloacae"
            "Enterococcus faecalis"
            "Haemophilus influenzae"    
            "Klebsiella oxytoca"
            "Neisseria meningitidis"
            "Neisseria gonorrhoeae"
            "Pseudomonas aeruginosa" 
            "Serratia marcescens"  
            "Staphylococcus aureus"
            "Staphylococcus pseudintermedius"
            "Streptococcus agalactiae"
            "Streptococcus pyogenes"
            "Vibrio cholerae"
            "Vibrio parahaemolyticus"
            "Vibrio vulnificus"
        )

        # check organism in curated organism list
        genus=$(echo "~{organism}" | cut -d " " -f1)
        taxon=$(echo "~{organism}" | cut -d " " -f1,2)
        
        amrfinder_organism=""
        if [[ "$genus" == "Acinetobacter" ]]; then
            for i in "${abcc[@]}"; do
                if [[ "$taxon" == "$i" ]]; then
                    amrfinder_organism="Acinetobacter_baumannii"
                    break
                fi
            done
        elif [[ "$genus" == "Burkholderia" ]]; then
            for i in "${bcc[@]}"; do
                if [[ "$taxon" == "$i" ]]; then
                    amrfinder_organism="Burkholderia_cepacia"
                break
                fi
            done
            for i in "${bpc[@]}"; do
                if [[ "$taxon" == "$i" ]]; then
                    amrfinder_organism="Burkholderia_pseudomallei"
                    break
                fi
            done
        elif [[ "$genus" == "Shigella" ]] || [[ "$genus" == "Escherichia" ]]; then
            amrfinder_organism="Escherichia"
        elif [[ "$genus" == "Salmonella" ]]; then
            amrfinder_organism="Salmonella"
        elif [[ "$taxon" == "Campylobacter coli" ]] || [[ "$taxon" == "Campylobacter jejuni" ]]; then
            amrfinder_organism="Campylobacter"
        elif [[ "$taxon" == "Enterococcus faecium" ]] || [[ "$taxon" == "Enterococcus hirae" ]]; then
            amrfinder_organism="Enterococcus_faecium"
        elif [[ "$taxon" == "Klebsiella pneumoniae" ]] || [[ "$taxon" == "Klebsiella aerogenes" ]]; then
            amrfinder_organism="Klebsiella_pneumoniae"
        elif [[ "$taxon" == "Streptococcus pneumoniae" ]] || [[ "$taxon" == "Streptococcus mitis" ]]; then
            amrfinder_organism="Streptococcus_pneumoniae"
        else    
            for i in "${taxa[@]}"; do
                if [[ "$taxon" == "$i" ]]; then
                    amrfinder_organism=${taxon// /_}
                    break
                fi
            done
        fi

        # checking bash variable
        echo "amrfinder organism is set to: ${amrfinder_organism}"

        # run amrfinderplus
        amrfinder \
        --threads ~{cpu} \
        --plus \
        --organism "$amrfinder_organism" \
        --name ~{id} \
        --nucleotide ~{assembly} \
        --protein ~{bakta_faa} \
        --gff ~{bakta_gff} \
        --annotation_format bakta \
        -o ~{id}.amrfinder.tsv
    >>>

    output {
        String amrfinder_version = read_string("VERSION")
        String amrfinder_db_version = read_string("DB_VERSION")
        File amrfinder_report = "~{id}.amrfinder.tsv"
    }

    runtime {
        docker: "staphb/ncbi-amrfinderplus:4.2.7-2026-03-24.1"
        cpu: cpu
        memory: "8GiB"
    }
}