version 1.0

import "../tasks/task_pbtk.wdl" as pbtk
import "../tasks/task_lrge.wdl" as lrge
import "../tasks/task_rasusa.wdl" as rasusa
import "../tasks/task_hifiasm.wdl" as hifiasm
import "../tasks/task_flye.wdl" as flye
import "../tasks/task_wtdbg2.wdl" as wtdbg2
import "../tasks/task_raven.wdl" as raven
import "../tasks/task_autocycler.wdl" as autocycler
import "../tasks/task_plassembler.wdl" as plassembler
import "../tasks/task_fastp.wdl" as fastp
import "../tasks/task_polypolish.wdl" as polypolish
import "../tasks/task_dnaapler.wdl" as dnaapler
import "../tasks/task_minimap2.wdl" as minimap2
import "../tasks/task_merge.wdl" as merge
import "../tasks/task_bandage.wdl" as bandage

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
    }

    input {
        String id
        File long_fq
        File? short_fq1
        File? short_fq2
    }

    # convert to fastq if bam file is provided
    String filename = basename(long_fq)

    if (sub(filename, ".bam", "") != filename) {
        call pbtk.bam_to_fastq {
            input:
                id = id,
                bam = long_fq
        }
    }

    # estimate genome size
    call lrge.estimate_genome_size {
        input: 
            long_fq = select_first([bam_to_fastq.long_fq,long_fq])
    }

    # downsample reads
    call rasusa.downsample {
        input:
            id = id,
            long_fq = select_first([bam_to_fastq.long_fq,long_fq]),
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

    # wtdbg2
    call wtdbg2.wtdbg2_asm {
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

    # autocycler
    call autocycler.combine_asms {
        input:
            id = id,
            flye_asm = flye_asm.assembly_fasta,
            hifiasm_asm = hifiasm_asm.assembly_fasta,
            wtdbg2_asm = wtdbg2_asm.assembly_fasta,
            raven_asm = raven_asm.assembly_fasta
    }

    # plasmid recovery & polishing
    if (defined(short_fq1) && defined(short_fq2)) {        
        # fastp
        call fastp.trim_pe {
            input:
                id = id,
                short_fq1 = short_fq1,
                short_fq2 = short_fq2
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
                autocycler_asm = combine_asms.assembly_fasta,
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
            long_asm = select_first([polish.polished_fasta, combine_asms.assembly_fasta])
    }

    call bandage.asm_image {
        input:
            id = id,
            hifiasm_gfa = hifiasm_asm.assembly_graph,
            flye_gfa = flye_asm.assembly_graph,
            raven_gfa = raven_asm.assembly_graph,
            wtdbg2_asm = wtdbg2_asm.assembly_fasta,
            autocycler_gfa = combine_asms.assembly_graph,
            plassembler_gfa = plassembler_asm.graph,
            final_asm = reorient.reoriented_fasta
    }

    # outputs
    output {
        # haphir version
        String version = "HAPHiR v0.2.0"
        # pbtk
        String? bam2fastq_version = bam_to_fastq.bam2fastq_version
        # lrge
        String lrge_version = estimate_genome_size.lrge_version
        String genome_size = estimate_genome_size.genome_size
        # flye
        String flye_version = flye_asm.flye_version
        File flye_assembly = flye_asm.assembly_fasta
        File flye_graph = flye_asm.assembly_graph
        # hifiasm
        String hifiasm_version = hifiasm_asm.hifiasm_version
        File hifiasm_assembly = hifiasm_asm.assembly_fasta
        File hifiasm_graph = hifiasm_asm.assembly_graph
        # wtdb2
        String wtdbg2_version = wtdbg2_asm.wtdbg2_version
        File wtdb2_assembly = wtdbg2_asm.assembly_fasta
        # raven
        String raven_version = raven_asm.raven_version
        File raven_assembly = raven_asm.assembly_fasta
        File raven_graph = raven_asm.assembly_graph
        # autocycler
        String autocycler_version = combine_asms.autocycler_version
        File autocycler_assembly = combine_asms.assembly_fasta
        File autocycler_graph = combine_asms.assembly_graph
        # fastp
        String? fastp_version = trim_pe.fastp_version
        File? fastp_report = trim_pe.html_report
        # plassembler
        String? plassembler_version = plassembler_asm.plassembler_version
        File? plassembler_plasmids = plassembler_asm.plasmids
        File? plassembler_graph = plassembler_asm.graph
        File? plassembler_summary = plassembler_asm.summary
        # minimap2
        String? minimap2_version = label_and_align.minimap_version
        File? minimap2_report = label_and_align.overlaps_paf
        # assembly merging 
        File? merged_fasta = merge_asms.merged_fasta
        File? merge_summary = merge_asms.merge_summary
        # polypolish
        String? polypolish_version = polish.polypolish_version
        File? polished_assembly = polish.polished_fasta
        # dnaapler
        String dnaapler_version = reorient.dnaapler_version
        File dnaapler_summary = reorient.dnaapler_summary
        File final_assembly = reorient.reoriented_fasta
        # bandage
        String bandage_version = asm_image.bandage_version
        File bandage_html = asm_image.bandage_html
    }
}