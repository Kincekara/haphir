version 1.1

import "wf_haphir.wdl" as main

workflow haphir_local {
    input {
        File samplesheet
        Array[Array[String]] inputSamples = read_tsv(samplesheet)        
    }

    scatter (sample in inputSamples) {
        if (length(sample) == 4) {
            call main.haphir as hybrid {            
                input:
                    id = sample[0],
                    long_fq = sample[1],
                    short_fq1 = sample[2],
                    short_fq2 = sample[3]
            }
        }
        if (length(sample) == 2) {
            call main.haphir as hifi {            
                input:
                    id = sample[0],
                    long_fq = sample[1]
            }                    
        }
    }
}
