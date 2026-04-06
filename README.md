# HAPHiR: Hybrid Assembly of PacBio HiFi and Illumina Reads

[![CI](https://github.com/yourusername/haphir/actions/workflows/check-wdl.yml/badge.svg)](https://github.com/yourusername/haphir/actions/workflows/check-wdl.yml)

***This repo is under development!***

HAPHiR performs high‑quality bacterial genome assembly using PacBio HiFi long reads and Illumina short reads, combining accuracy, robustness, and efficient cloud execution.

The workflow runs multiple long‑read assemblers in parallel (Flye, Hifiasm, Raven, wtdbg2) and generates a unified, high‑confidence consensus assembly using [Autocycler](https://github.com/rrwick/Autocycler). Small circular plasmids are recovered through a dedicated hybrid assembly step using [Plassembler](https://github.com/gbouras13/plassembler), ensuring both chromosomal and plasmid components are accurately reconstructed.

HAPHiR is designed for cloud‑native execution on [Terra](https://terra.bio/), but can also be run locally using WDL executer such as [miniwdl](https://miniwdl.readthedocs.io/en/latest/) or [Cromwell](https://cromwell.readthedocs.io/en/latest/).

## Features

- **Multi-assembler consensus**: Runs 4 independent long-read assemblers (Flye, Hifiasm, Wtdbg2, Raven) and combines them using Autocycler for enhanced accuracy
- **HiFI only support**: Works with PacBio HiFi-only data or hybrid HiFi + Illumina data
- **Plasmid recovery**: Dedicated plasmid assembly and recovery using Plassembler
- **Flexible inputs**: Accepts PacBio BAM or FASTQ files, automatically converts as needed
- **Quality control**: Includes read trimming, genome size estimation, and coverage normalization
- **Polishing**: Optional short-read polishing with Polypolish 
- **Cloud-ready**: Designed for scalable execution on Terra
- **Containerized**: All tools run in Docker containers for reproducibility

## Quick Start

### Prerequisites


### Single Sample Assembly


### Batch Processing

## Input Files

### Single Sample Workflow (`wf_haphir.wdl`)

| Input | Type | Description |
|-------|------|-------------|
| `id` | String | Sample identifier |
| `long_fq` | File | PacBio HiFi reads (FASTQ or BAM) |
| `short_fq1` | File (optional) | Illumina forward reads |
| `short_fq2` | File (optional) | Illumina reverse reads |
| `polishing` | Boolean | Enable Illumina polishing (default: true) |

### Batch Workflow (`wf_haphir_local.wdl`)

Provide a TSV file with sample information:

```
# HiFi-only (2 columns)
sample1	path/to/sample1.hifi.fastq
sample2	path/to/sample2.hifi.fastq

# Hybrid (4 columns)
sample3	path/to/sample3.hifi.fastq	path/to/sample3.R1.fastq	path/to/sample3.R2.fastq
```

## Pipeline Overview

1. **Input Processing**: Convert BAM to FASTQ if needed
2. **Genome Characterization**: Estimate genome size and normalize coverage to  ~100X
3. **Parallel Assembly**: Run 4 long-read assemblers simultaneously
4. **Consensus Generation**: Combine assemblies using Autocycler
5. **Reorientation**: Orient assembly to DnaA origin using Dnaapler
6. **Polishing**: Optional short-read polishing with Polypolish
7. **Plasmid Recovery**: Assemble plasmids using Plassembler
8. **Final Assembly**: Merge chromosome and plasmids

## Output Files

All outputs are prefixed with the sample ID:

- `{id}.final_assembly.fasta` - Complete assembly (chromosome + plasmids)
- `{id}.seqkit.stats.txt` - Assembly statistics
- `{id}.plasmids.fasta` - Recovered plasmid sequences
- `{id}.polished.fasta` - Polished assembly (if polishing enabled)
- `{id}.dnaapler.fasta` - Reoriented assembly
- `{id}.autocycler.fasta` - Consensus assembly
- Individual assembler outputs and intermediate files

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Flye | 2.9.6 | Long-read assembler |
| HiFiASM | 0.25.0 | HiFi-optimized assembler |
| Wtdbg2 | 2.5 | De Bruijn graph assembler |
| Raven | 1.8.3 | Modern graph-based assembler |
| Autocycler | 0.6.1 | Assembly consensus |
| Dnaapler | 1.3.0 | Assembly reorientation |
| Polypolish | 0.6.1 | Short-read polishing |
| Plassembler | 1.8.2 | Plasmid assembly |
| Fastp | 1.1.0 | Read trimming |
| Seqkit | 2.13.0 | Sequence manipulation |
| LRGE | 0.2.1 | Genome size estimation |
| Rasusa | 3.0.0 | Coverage normalization |

## Computing Requirements

- **Hifiasm**: 16 CPU, 32GB RAM (most resource-intensive)
- **Other assemblers**: 8 CPU, 16GB RAM
- **Storage**: 200GB SSD recommended for large genomes
- **Runtime**: 2-8 hours depending on genome size and compute resources

## Configuration

The pipeline is designed to run on various platforms:

- **Local execution**: Miniwdl + Docker
- **Cloud platforms**: Terra
- **HPC clusters**: Cromwell + Singularity

## Validation

The repository includes automated validation:

- WDL syntax checking with MiniWDL
- Cromwell workflow validation
- CI/CD pipeline for pull requests

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add/update tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use HAPHiR in your research, please cite:

```

```

## Support

For issues, questions, or feature requests, please open an issue on GitHub.