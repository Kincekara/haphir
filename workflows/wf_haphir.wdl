version 1.0

import "../tasks/task_lrge.wdl" as lrge
import "../tasks/task_rasusa.wdl" as rasusa
import "../tasks/task_hifiasm.wdl" as hifiasm
import "../tasks/task_flye.wdl" as flye
import "../tasks/task_wtdbg2.wdl" as wtdbg2
import "../tasks/task_raven.wdl" as raven
import "../tasks/task_canu.wdl" as canu
import "../tasks/task_lja.wdl" as lja
import "../tasks/task_autocycler.wdl" as autocycler
import "../tasks/task_plassembler.wdl" as plassembler
import "../tasks/task_fastp.wdl" as fastp
import "../tasks/task_polypolish.wdl" as polypolish
import "../tasks/task_dnaapler.wdl" as dnaapler
import "../tasks/task_minimap2.wdl" as minimap2
import "../tasks/task_merge.wdl" as merge
import "../tasks/task_bandage.wdl" as bandage
import "../tasks/task_bakta.wdl" as bakta
import "../tasks/task_amrfinderplus.wdl" as amrfinderplus

workflow haphir {
    meta {
        author: "Kutluhan Incekara"
        email: "kutluhan.incekara@ct.gov"
        description: "Hybrid Assembly of PacBio HiFi and Illumina Reads"
    }

    parameter_meta {
        id: {
            description: "Unique identifier for the assembly"
        }
        long_fq: {
            description: "PacBio HiFi reads in bam or fastq format",
            patterns: [".bam", ".fastq.gz"]
        }
        short_fq1: {
            description: "Illumina paired-end reads 1",
            patterns: ["*_1.fastq.gz", "*_R1_*.fastq.gz"]
        }
        short_fq2: {
            description: "Illumina paired-end reads 2",
            patterns: ["*_2.fastq.gz", "*_R2_*.fastq.gz"]
        }
        organism: {
            description: "taxonomic name of the organism",
            patterns: ["^[A-Z][a-z]+ [a-z]+$", "^[A-Z][a-z]+ [a-z]+( [a-z]+\.? [a-z0-9-]+)?$"]
        }
        bakta_annotation: {
            description: "Run bakta for annotation",
            patterns: ["true", "false"],
            default: "false"
        }
        amrfinder: {
            description: "Run amrfinder for antibiotic resistance gene detection",
            patterns: ["true", "false"],
            default: "false"
        }
        assembly_preset: {
            description: "Assembly preset to determine the order of assembly tools",
            patterns: ["1", "2", "3"],
            default: "1"
        }
    }

    input {
        String id
        File long_fq
        File? short_fq1
        File? short_fq2
        String? organism
        Boolean bakta_annotation = false
        Boolean amrfinder = false
        Int assembly_preset = 1
    }

    # Presets    
    Array[Array[String]] assembly_presets = [["wtdbg2", "canu", "lja"],["canu", "lja", "wtdbg2"],["lja", "wtdbg2", "canu"]]
    Array[String] assembly_option = assembly_presets[assembly_preset - 1]
    
    # estimate genome size
    call lrge.estimate_genome_size {
        input: 
            long_fq = long_fq
    }

    # downsample reads
    call rasusa.downsample {
        input:
            id = id,
            long_fq = long_fq,
            genome_size = estimate_genome_size.genome_size
    }

    # flye
    call flye.flye_asm {
        input:
            id = id,
            long_fq = downsample.downsampled_fq,
            genome_size = estimate_genome_size.genome_size
    }

    # hifiasm
    call hifiasm.hifiasm_asm {
        input:
            id = id,
            long_fq = downsample.downsampled_fq,
            genome_size = estimate_genome_size.genome_size
    }

    # raven
    call raven.raven_asm {
        input:
            id = id,
            long_fq = downsample.downsampled_fq
    }
    
    # ==================== OPTION 1 / ATTEMPT 1 ====================
    # wtdbg2
    if (assembly_option[0] == "wtdbg2") {
        call wtdbg2.wtdbg2_asm as wtdbg2_opt1 {
            input:
                id = id,
                long_fq = downsample.downsampled_fq,
                genome_size = estimate_genome_size.genome_size
        }
    }

    if (assembly_option[0] == "canu") {
        call canu.hicanu_asm as hicanu_opt1 {
            input:
                id = id,
                long_fq = downsample.downsampled_fq,
                genome_size = estimate_genome_size.genome_size
        }
    }
    
    if (assembly_option[0] == "lja") {
        call lja.lja_asm as lja_opt1 {
            input:
                id = id,
                long_fq = downsample.downsampled_fq
        }
    }

    File opt1_asm = select_first([wtdbg2_opt1.assembly_fasta, hicanu_opt1.assembly_fasta, lja_opt1.assembly_fasta])

    # autocycler
    call autocycler.combine_asms as autocycler_attempt1 {
        input:
            id = id,
            flye_asm = flye_asm.assembly_fasta,
            hifiasm_asm = hifiasm_asm.assembly_fasta,            
            raven_asm = raven_asm.assembly_fasta,
            opt_asm = opt1_asm
    }

    # ==================== OPTION 2 / ATTEMPT 2 ====================
    if (autocycler_attempt1.autocycler_result == "FAIL") {
        # wtdbg2
        if (assembly_option[1] == "wtdbg2") {
            call wtdbg2.wtdbg2_asm as wtdbg2_opt2 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq,
                    genome_size = estimate_genome_size.genome_size
            }
        }

        if (assembly_option[1] == "canu") {
            call canu.hicanu_asm as hicanu_opt2 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq,
                    genome_size = estimate_genome_size.genome_size
            }
        }
    
        if (assembly_option[1] == "lja") {
            call lja.lja_asm as lja_opt2 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq
            }
        }

        File opt2_asm = select_first([wtdbg2_opt2.assembly_fasta, hicanu_opt2.assembly_fasta, lja_opt2.assembly_fasta])

        # autocycler
        call autocycler.combine_asms as autocycler_attempt2 {
            input:
                id = id,
                flye_asm = flye_asm.assembly_fasta,
                hifiasm_asm = hifiasm_asm.assembly_fasta,            
                raven_asm = raven_asm.assembly_fasta,
                opt_asm = opt2_asm
        }
    }
    # ==================== OPTION 3 / ATTEMPT 3 ====================
    if ((autocycler_attempt1.autocycler_result == "FAIL") && (autocycler_attempt2.autocycler_result == "FAIL")) {
        # wtdbg2
        if (assembly_option[2] == "wtdbg2") {
            call wtdbg2.wtdbg2_asm as wtdbg2_opt3 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq,
                    genome_size = estimate_genome_size.genome_size
            }
        }

        if (assembly_option[2] == "canu") {
            call canu.hicanu_asm as hicanu_opt3 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq,
                    genome_size = estimate_genome_size.genome_size
            }
        }
    
        if (assembly_option[2] == "lja") {
            call lja.lja_asm as lja_opt3 {
                input:
                    id = id,
                    long_fq = downsample.downsampled_fq
            }
        }

        File opt3_asm = select_first([wtdbg2_opt3.assembly_fasta, hicanu_opt3.assembly_fasta, lja_opt3.assembly_fasta])

        # autocycler
        call autocycler.combine_asms as autocycler_attempt3 {
            input:
                id = id,
                flye_asm = flye_asm.assembly_fasta,
                hifiasm_asm = hifiasm_asm.assembly_fasta,            
                raven_asm = raven_asm.assembly_fasta,
                opt_asm = opt3_asm
        }
    }
    # ==================== END OF AUTOCYCLER ATTEMPTS ====================   
    
    # prep short reads
    if (defined(short_fq1) && defined(short_fq2)) { 
        # rasusa
        call rasusa.downsample_pe {
            input:
                id = id,
                short_fq1 = short_fq1,
                short_fq2 = short_fq2,
                genome_size = estimate_genome_size.genome_size
        }
        # fastp
        call fastp.trim_pe {
            input:
                id = id,
                short_fq1 = downsample_pe.ds_short_fq1,
                short_fq2 = downsample_pe.ds_short_fq2
        }
    }

    # plassembler
    call plassembler.plassembler_asm {
        input:
            id = id,
            long_fq = downsample.downsampled_fq,
            short_fq1 = trim_pe.short_fq1_trimmed,
            short_fq2 = trim_pe.short_fq2_trimmed,
            flye_asm = flye_asm.assembly_fasta,
            flye_info = flye_asm.assembly_info
    }
    
    # minimap2
    call minimap2.label_and_align {
        input:
            id = id,
            autocycler_asm = select_first([autocycler_attempt3.assembly_fasta, autocycler_attempt2.assembly_fasta, autocycler_attempt1.assembly_fasta]),
            plassembler_asm = plassembler_asm.plasmids
    }

    # merge assemblies
    call merge.merge_asms {
        input:
            id = id,
            autocycler_asm = label_and_align.autocycler_fasta,
            plassembler_asm = label_and_align.plasmids_fasta,
            overlaps_paf = label_and_align.overlaps_paf
    }

    # polypolish
    if (defined(short_fq1) && defined(short_fq2)) { 
        call polypolish.polish {
            input:
                id = id,
                draft_asm = merge_asms.merged_fasta,
                short_fq1 = trim_pe.short_fq1_trimmed,
                short_fq2 = trim_pe.short_fq2_trimmed
        }
    }

    # dnaapler
    call dnaapler.reorient {
        input:
            id = id,
            long_asm = select_first([polish.polished_fasta, merge_asms.merged_fasta])
    }

    # bandage
    call bandage.asm_image {
        input:
            id = id,
            hifiasm_gfa = hifiasm_asm.assembly_graph,
            flye_gfa = flye_asm.assembly_graph,
            raven_gfa = raven_asm.assembly_graph,
            opt_gfa = select_first([wtdbg2_opt3.assembly_fasta, hicanu_opt3.assembly_fasta, lja_opt3.assembly_graph,
                                    wtdbg2_opt2.assembly_fasta, hicanu_opt2.assembly_fasta, lja_opt2.assembly_graph,
                                    wtdbg2_opt1.assembly_fasta, hicanu_opt1.assembly_fasta, lja_opt1.assembly_graph]),
            autocycler_gfa = select_first([autocycler_attempt3.assembly_graph, autocycler_attempt2.assembly_graph, autocycler_attempt1.assembly_graph]),
            plassembler_gfa = plassembler_asm.graph,
            final_asm = reorient.reoriented_fasta,
            hifiasm_ctg_len = hifiasm_asm.ctg_len,
            flye_ctg_len = flye_asm.ctg_len,
            raven_ctg_len = raven_asm.ctg_len,
            opt_ctg_len = select_first([wtdbg2_opt3.ctg_len, hicanu_opt3.ctg_len, lja_opt3.ctg_len,
                                        wtdbg2_opt2.ctg_len, hicanu_opt2.ctg_len, lja_opt2.ctg_len,
                                        wtdbg2_opt1.ctg_len, hicanu_opt1.ctg_len, lja_opt1.ctg_len]),
            autocycler_ctg_len = select_first([autocycler_attempt3.ctg_len, autocycler_attempt2.ctg_len, autocycler_attempt1.ctg_len]),
            plassembler_ctg_len = plassembler_asm.ctg_len,
            final_ctg_len = reorient.ctg_len
        }

    if ( bakta_annotation || amrfinder ) {
        call bakta.annotation {
            input:
                id = id,
                final_asm = reorient.reoriented_fasta,
                organism = organism
        }
    }

    if ( amrfinder ) {
        call amrfinderplus.amr {
            input:
                id = id,
                assembly = annotation.bakta_fna,
                bakta_faa = annotation.bakta_faa,
                bakta_gff = annotation.bakta_gff,
                organism = organism
        }
    }

    # outputs
    output {
        # haphir version
        String version = "HAPHiR v0.12.0"
        # rasusa
        String est_longfq_cov = downsample.coverage
        String? est_shortfq_cov = downsample_pe.coverage
        # autocycler
        File autocycler_assembly = select_first([autocycler_attempt3.assembly_fasta, autocycler_attempt2.assembly_fasta, autocycler_attempt1.assembly_fasta])
        File autocycler_graph = select_first([autocycler_attempt3.assembly_graph, autocycler_attempt2.assembly_graph, autocycler_attempt1.assembly_graph])
        # fastp
        File? fastp_report = trim_pe.html_report
        # plassembler
        File plassembler_plasmids = plassembler_asm.plasmids
        File plassembler_graph = plassembler_asm.graph
        File plassembler_summary = plassembler_asm.summary
        # minimap2
        File minimap2_report = label_and_align.overlaps_paf
        # assembly merging 
        File merge_summary = merge_asms.merge_summary
        # dnaapler
        File dnaapler_summary = reorient.dnaapler_summary
        File final_assembly = reorient.reoriented_fasta
        # bandage
        File asm_viz = asm_image.bandage_html
        # bakta
        File? bakta_outputs = annotation.bakta_outputs
        # amrfinderplus
        File? amrfinder_report = amr.amrfinder_report
        # program versions
        Array[String] program_versions = [ "lrge: " + estimate_genome_size.lrge_version,
                                "rasusa: " + downsample.rasusa_version,
                                "flye: " + flye_asm.flye_version,
                                "hifiasm: " + hifiasm_asm.hifiasm_version,
                                "wtdbg2: " + select_first([wtdbg2_opt3.wtdbg2_version, wtdbg2_opt2.wtdbg2_version, wtdbg2_opt1.wtdbg2_version, "NA"]),
                                "canu: " + select_first([hicanu_opt3.canu_version, hicanu_opt2.canu_version, hicanu_opt1.canu_version, "NA"]),
                                "lja: " + select_first([lja_opt3.lja_version, lja_opt2.lja_version, lja_opt1.lja_version, "NA"]),
                                "raven: " + raven_asm.raven_version,
                                "autocycler: " + select_first([autocycler_attempt3.autocycler_version, autocycler_attempt2.autocycler_version, autocycler_attempt1.autocycler_version, "NA"]),
                                "fastp: " + select_first([trim_pe.fastp_version, "NA"]),
                                "plassembler: " + plassembler_asm.plassembler_version,
                                "minimap2: " + label_and_align.minimap_version,
                                "polypolish: " + select_first([polish.polypolish_version, "NA"]),
                                "dnaapler: " + reorient.dnaapler_version,
                                "bandage: " + asm_image.bandage_version,
                                "bakta: " + select_first([annotation.bakta_version, "NA"]),
                                "amrfinderplus: " + select_first([amr.amrfinder_version, "NA"])
                                ]
        
    }
}