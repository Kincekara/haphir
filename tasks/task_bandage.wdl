version 1.0

task asm_image {
    input {
        String id
        File hifiasm_gfa
        File flye_gfa
        File raven_gfa
        File wtdbg2_asm
        File autocycler_gfa
        File? plassembler_gfa
        File final_asm
    }

    command <<< 
        set -euo pipefail

        # version
        Bandage --version | cut -d " " -f2 > VERSION

        # Bandage
        Bandage image ~{hifiasm_gfa} hifiasm.png
        Bandage image ~{flye_gfa} flye.png
        Bandage image ~{raven_gfa} raven.png
        Bandage image ~{wtdbg2_asm} wtdbg2.png
        Bandage image ~{autocycler_gfa} autocycler.png
        if [ -n "~{plassembler_gfa}" ]; then
            Bandage image ~{plassembler_gfa} plassembler.png
        fi
        Bandage image ~{final_asm} final.png

        # write html file
        cat << EOF > ~{id}.bandage.html
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
                </style>            
            </head>
            <body>
                <h1>Assembly Comparison</h1>
                <h2>Intermediate Assemblies</h2>
                <div class="grid">
                    <div class="grid-item">
                        <div class="caption">Hifiasm</div>
                        <img src="data:image/png;base64,$(base64 -w 0 hifiasm.png)" alt="Hifiasm" />
                    </div>
                    <div class="grid-item">
                        <div class="caption">Flye</div>
                        <img src="data:image/png;base64,$(base64 -w 0 flye.png)" alt="Flye" />
                    </div>
                    <div class="grid-item">
                        <div class="caption">Raven</div>
                        <img src="data:image/png;base64,$(base64 -w 0 raven.png)" alt="Raven" />
                    </div>
                    <div class="grid-item">
                        <div class="caption">Wtdbg2</div>
                        <img src="data:image/png;base64,$(base64 -w 0 wtdbg2.png)" alt="Wtdbg2" />
                    </div>
                    <div class="grid-item">
                        <div class="caption">Autocycler</div>
                        <img src="data:image/png;base64,$(base64 -w 0 autocycler.png)" alt="Autocycler" />
                    </div>
                    <!--<div class="grid-item">
                        <div class="caption">Plassembler</div>
                        <img src="data:image/png;base64,$(base64 -w 0 plassembler.png)" alt="Plassembler">
                    </div>-->
                </div>
                <br><br>
                <h2>Final Assembly</h2>
                <img src="data:image/png;base64,$(base64 -w 0 final.png)" alt="final_asm">                
            </body>
            <hr>
            <footer>
            <p><i>This report is created by <a href="https://github.com/Kincekara/haphir">HAPHiR</a> bioinformatics pipeline.</br></i></p>
            </footer>
        </html>
        EOF

        if [ -f plassembler.png ]; then
            sed -i 's/<!--//g; s/-->//g' ~{id}.bandage.html
        fi
    >>>

    output {
        String bandage_version = read_string("VERSION")
        File bandage_html = "~{id}.bandage.html"
    }

    runtime {
        docker: "kincekara/bandage:0.9.0"
        cpu: 1
        memory: "1 GiB"
    }
}