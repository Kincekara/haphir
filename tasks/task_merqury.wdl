version 1.0

task asm_qc {
    input {
        String id
        File asm
        File? long_fq
        File? short_fq1
        File? short_fq2
        Int cpu = 4
    }

    command <<<
        set -euo pipefail        

        # meryl - use short reads if available, otherwise use long reads
        if [[ -n "~{short_fq1}" && -n "~{short_fq2}" ]]; then
            mode="hybrid"
            meryl count threads=~{cpu} k=19 ~{short_fq1} ~{short_fq2} output reads.meryl
        else
            mode="hifi"
            meryl count threads=~{cpu} k=19 ~{long_fq} output reads.meryl
        fi   

        # merqury
        merqury.sh reads.meryl ~{asm} ~{id}

        # parse results
        awk -F'\t' '{print $5}' ~{id}.completeness.stats > COMPLETENESS
        awk -F'\t' '{print $4}' ~{id}.qv > QV

        # evaluate assembly quality
        comp=$(< COMPLETENESS)
        qv=$(< QV)

        if [[ "$mode" == "hybrid" ]]; then
            if awk -v comp="$comp" -v qv="$qv" \
                'BEGIN {
                    if (comp > 99.0 && (qv == "inf" || qv > 40.0)) exit 0
                    exit 1
                }'
            then
                echo "PASS" > EVAL
            else
                echo "WARN[C:$comp,QV:$qv]" > EVAL
            fi
        fi

        if [[ "$mode" == "hifi" ]]; then
            if awk -v comp="$comp" -v qv="$qv" \
                'BEGIN {
                    if (comp > 98.5 && (qv == "inf" || qv > 35.0)) exit 0
                    exit 1
                }'
            then
                echo "PASS" > EVAL
        else
                echo "WARN[C:$comp,QV:$qv]" > EVAL
            fi
        fi

        # create table of completeness and qv
        awk -F'\t' '
            BEGIN { 
                print "  <tr>" 
            }
            {
                for (i = 1; i <= NF; i++) {
                    print "    <td>" $i "</td>"
                }
            }
            END { 
                print "  </tr>" 
            }' ~{id}.qv > qv_table

        awk -F'\t' '
            BEGIN { 
                print "  <tr>" 
            }
            {
                for (i = 1; i <= NF; i++) {
                    print "    <td>" $i "</td>"
                }
            }
            END { 
                print "  </tr>" 
            }' ~{id}.completeness.stats > comp_table


        # create html report
        cat << EOF > ~{id}.merqury.html
        <html>
            <head>
                <title>~{id}</title>
                <meta charset="utf-8" />
                <style>
                    body { font-family: Arial, sans-serif; margin: 20px; }
                    h1 { margin-bottom: 16px; text-align: center; }
                    .grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px; width: 100%; }
                    .grid-item { border: 1px solid #ccc; padding: 10px; }
                    .grid-item .caption { margin-top: 8px; font-weight: bold; text-align: center; }
                    img { max-width: 100%; height: auto; display: block; margin: 0 auto; }
                    table {font-family: arial, sans-serif; border-collapse: collapse; }
                    td, th { border: 1px solid #dddddd; text-align: left; padding: 6px; }
                    tr:nth-child(even) {  background-color: #dddddd; }
                </style>            
            </head>
            <body>
                <h1>Assembly Quality Control Report</h1>
                <h2>Reference free QV estimation</h2>
                <table>
                    <tr>
                        <th>Assembly</th>
                        <th>Num. of unique k-mers</th>
                        <th>Num. of k-mers</th>
                        <th>QV</th>
                        <th>Error rate</th>
                    </tr>
                    $(cat qv_table)
                </table>
                <br><br>

                <h2>K-mer Completeness</h2>
                <table>
                    <tr>
                        <th>Assembly</th>
                        <th>Read set</th>
                        <th>Num. of solid k-mers (asm)</th>
                        <th>Num. of solid k-mers (reads)</th>
                        <th>Completeness (%)</th>
                    </tr>
                    $(cat comp_table)
                </table>
                <br><br>
                
                <h2>Copy-number spectrum</h2>
                <div class="grid">
                    <div class="grid-item">
                        <div class="caption">stacked</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.~{id}.spectra-cn.st.png)" alt=stacked />
                    </div>
                    <div class="grid-item">
                        <div class="caption">filled</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.~{id}.spectra-cn.fl.png)" alt=filled />
                    </div>
                    <div class="grid-item">
                        <div class="caption">line</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.~{id}.spectra-cn.ln.png)" alt=line />
                    </div>
                </div>
                <h2>Assembly spectrum</h2>
                <div class="grid">
                    <div class="grid-item">
                        <div class="caption">stacked</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.spectra-asm.st.png)" alt=stacked />
                    </div>
                    <div class="grid-item">
                        <div class="caption">filled</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.spectra-asm.fl.png)" alt=filled />
                    </div>
                    <div class="grid-item">
                        <div class="caption">line</div>
                        <img src="data:image/png;base64,$(base64 -w 0 ~{id}.spectra-asm.ln.png)" alt=line />
                    </div>
                </div>
                <br><br>
            </body>
            <hr>
            <footer>
            <p><i>This report is created by <a href="https://github.com/Kincekara/haphir">HAPHiR</a> bioinformatics pipeline.</br></i></p>
            </footer>
        </html>
        EOF
    >>>

    output {
        String merqury_version = "1.4.1"
        File merqury_log = "logs/~{id}.spectra-cn.log"
        String eval = read_string("EVAL")
        String completeness = read_string("COMPLETENESS")
        String qv = read_string("QV")
        File merqury_html = "~{id}.merqury.html"

    }

    runtime {
        docker: "kincekara/merqury:1.4.1"
        cpu: cpu
        memory: "8 GiB"
        preemptible: 0
        maxRetries: 3
    }
}
