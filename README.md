# metagenomicsPipeline
A pipeline for MAG construction and species abundance estimation from metagenomic data

## Workflow
- The pipeline starts with raw shotgun sequencing reads, followed by trimming and quality filtering.
- Rough species abundance estimation with MetaPhlan and pathway abundance estimation with HUMANN
- Assembly by Spades and Megahit in parallel
- Binning with BASALT
- Bin dereplication using dRep
- mapping of reads against the non-redundant bin catalogue with Bowtie2 and read counting with samtools
- bin annotation with Prokka
- construction of a protein catalogue from the annotations, followed by 100% identity clustering and annotation via eggNOG
- Taxonomy assignment with gtdbTk

## Installation and notes
The pipeline is written in Nextflow and dependencies are managed with conda. We are currently not using the native conda interface from Nextflow yet, but rather activate environments manually. Environment files can be found in the respective folder.
The pipeline is split in two parts, one taking care of all steps until binning, while the second part is responsbile for the rest.

The second pipeline takes in a metaData.csv, which contains pairs of samples and individuals from which the samples originated. This can be usefull when not wanting to dereplicte across individuals to preserve potential strain differneces. The format should look like:

```
S1,IDV1
S2,IDV1
S3,IDV2
S4,IDV2
S5,IDV2
```

To run the pipeline:
```
nextflow ./01.CleanProfileBin.nf --reads './FASTQ/P137_*_{1,2}.fq.gz'
nextflow ./02.ClassAnnotMap.nf --basaltBins 'BASALT_BINS/*' --metaData metaData.csv --cleanReads 'CLEANREADS/P137_*_{1,2}.clean.fq.gz'
```
