version 1.0

import "wf_haphir.wdl" as main

workflow batch {
    input {
        File samplesheet
        Array[Array[String]] inputSamples = read_tsv(samplesheet)
        Boolean bakta = false
        Boolean amrfinder = false      
    }

    scatter (idx in range(length(inputSamples))) {

        if ( idx > 0 ) { #skip header

            #define columns
            String id = inputSamples[idx][0]
            String long_fq = inputSamples[idx][1]
            String short_fq1 = inputSamples[idx][2]
            String short_fq2 = inputSamples[idx][3]
            String organism = inputSamples[idx][4]

            if ( short_fq1 != "" && short_fq2 != "" ) {
                call main.haphir as hybrid {            
                    input:
                        id = id,
                        long_fq = long_fq,
                        short_fq1 = short_fq1,
                        short_fq2 = short_fq2,
                        organism = organism,
                        bakta_annotation = bakta,
                        amrfinder = amrfinder
                }
            }

            if ( short_fq1 == "" || short_fq2 == "" ) {
                call main.haphir as hifi {            
                    input:
                        id = id,
                        long_fq = long_fq,
                        organism = organism,
                        bakta_annotation = bakta,
                        amrfinder = amrfinder
                }  
            }
        }
    }
}
