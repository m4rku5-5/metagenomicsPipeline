#!/bin/nextflow

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
reads.view()

process PREP {
    cpus 24
    time '4h'
    memory '42 GB'
    scratch true
    stageOutMode 'move'

    publishDir "./CLEANREADS", mode: 'symlink'
    
    input:
         tuple val(sample_id), path(reads)

    output:
         tuple val(sample_id), path("${sample_id}_{1,2}.clean.fq.gz"), emit: cleanReads
    
    script:
    """
    conda activate PREPQC

    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_1.trim.fq.gz -O ${sample_id}_2.trim.fq.gz --dedup --dup_calc_accuracy 6 -l 50 -w 18

    hocort map biobloom -x DB/GRCh38/reference.bf -i ${sample_id}_1.trim.fq.gz ${sample_id}_2.trim.fq.gz -o ${sample_id}.HR -t 24

    mv ${sample_id}.HR_noMatch_1.fq ${sample_id}_1.clean.fq
    mv ${sample_id}.HR_noMatch_2.fq ${sample_id}_2.clean.fq
    pigz ${sample_id}_*.clean.fq
    """
}

process PROFILE_HUMANN {
    cpus 32
    time { 20.hour * task.attempt }
    memory { 96.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 2
    scratch true
    stageOutMode 'move'

    publishDir "./PROFILE/HUMANN", mode: 'copy', overwrite: true
    
    input:
         tuple val(sample_id), path(cleanReads)

    output:
         path("${sample_id}.comb_*.tsv")
    
    script:

    """
    conda activate PROFILE

    cat ${cleanReads[0]} ${cleanReads[1]} > ${sample_id}.comb.fq.gz
    humann --input ${sample_id}.comb.fq.gz --output ${sample_id} --threads 32
    
    mv ${sample_id}/${sample_id}.comb_genefamilies.tsv .
    mv ${sample_id}/${sample_id}.comb_pathabundance.tsv .
    mv ${sample_id}/${sample_id}.comb_pathcoverage.tsv .
    """
}

process PROFILE_METAPHLAN {
    cpus 24
    time { 4.hour * task.attempt }
    memory '48 GB'
    errorStrategy 'retry'
    maxRetries 2
    scratch true
    stageOutMode 'move'

    publishDir "./PROFILE/METAPHLAN", mode: 'copy', overwrite: true
    
    input:
         tuple val(sample_id), path(cleanReads)

    output:
         path("${sample_id}.metaphlan.tsv")
    
    script:

    conda activate PROFILE

    cat ${cleanReads[0]} ${cleanReads[1]} > ${sample_id}.comb.fq.gz
    metaphlan --input_type fastq ${sample_id}.comb.fq.gz -o "${sample_id}".metaphlan.tsv --nproc 24
    rm ${sample_id}.comb.fq.gz
    """
}

process PROFILE_MERGE_METAPHLAN {
    cpus 2
    time '10min'
    memory '4 GB'
    errorStrategy 'ignore'

    publishDir "./PROFILE/", mode: 'copy', overwrite: true
    
    input:
         path("*")

    output:
         path("merged.metaphlan.tsv")
    
    script:

    """
    conda activate PROFILE

    merge_metaphlan_tables.py -o merged.metaphlan.tsv *.metaphlan.tsv

    """
}

process PROFILE_MERGE_HUMANN {
    cpus 4
    time '30min'
    memory '24 GB'
    errorStrategy 'ignore'

    publishDir "./PROFILE/", mode: 'copy', overwrite: true
    
    input:
         path("*")

    output:
         path("merged_genefamilies.humann.tsv")
         path("merged_pathabundance.humann.tsv")
         path("merged_pathcoverage.humann.tsv")
         path("merged_genefamilies-cpm.humann.tsv")
         path("merged_pathabundance-cpm.humann.tsv")
         path("merged_pathcoverage-cpm.humann.tsv")
    
    script:

    """
    conda activate PROFILE

    humann_join_tables -i . -o merged_genefamilies.humann.tsv --file_name genefamilies
    humann_join_tables -i . -o merged_pathabundance.humann.tsv --file_name pathabundance
    humann_join_tables -i . -o merged_pathcoverage.humann.tsv --file_name pathcoverage
    
    humann_renorm_table -i merged_genefamilies.humann.tsv  -o merged_genefamilies-cpm.humann.tsv --units cpm
    humann_renorm_table -i merged_pathabundance.humann.tsv -o merged_pathabundance-cpm.humann.tsv --units cpm
    humann_renorm_table -i merged_pathcoverage.humann.tsv  -o merged_pathcoverage-cpm.humann.tsv --units cpm

    """
}

process ASM_MEGAHIT {
    cpus 24
    time { 5.h * task.attempt }
    memory { 64.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 2
    scratch true
    stageOutMode 'move'

    publishDir "./ASSEMBLY", mode: 'copy', overwrite: true
    
    input:
         tuple val(sample_id), path(cleanReads)

    output:
         path("${sample_id}.contigs.MEGAHIT.fa")
    
    script:

    """
    conda activate ASM

    megahit -1 ${cleanReads[0]} -2 ${cleanReads[1]} -t 24 -m 64e9 -o output

    mv output/final.contigs.fa ${sample_id}.contigs.MEGAHIT.fa
    """
}

process ASM_SPADES {
    cpus 32
    time { 20.h * task.attempt }
    memory { 196.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 2
    scratch true
    stageOutMode 'move'

    publishDir "./ASSEMBLY", mode: 'copy', overwrite: true
    
    input:
         tuple val(sample_id), path(cleanReads)

    output:
         path("${sample_id}.contigs.SPADES.fa")
    
    script:

    """
    conda activate ASM
	
    spades.py -1 ${cleanReads[0]} -2 ${cleanReads[1]} --meta -o output -t 32 -m 128

    mv output/scaffolds.fasta ${sample_id}.contigs.SPADES.fa
    """
}

process BASALT {
    cpus 32
    time { 32.h * task.attempt }
    memory 64.GB
    scratch true
    stageOutMode 'move'
    errorStrategy 'retry'
    maxRetries 2


    publishDir "./BASALT_BINS", mode: 'copy', overwrite: true
    
    input:
         tuple val(sample_id), path(cleanReads)
         path(assembly_megahit)
         path(assembly_spades)

    output:
         tuple val(sample_id), path("${sample_id}.FinalBinSet.tar.gz"), emit: basalt_bins
    
    script:

    """
    conda activate BASALT
      
    cp ${assembly_megahit} ./1_${sample_id}_megahit.fa
    cp ${assembly_spades} ./2_${sample_id}_spades.fa
    
    cp ${cleanReads[0]} ./${sample_id}_1.fq.gz
    cp ${cleanReads[1]} ./${sample_id}_2.fq.gz
    
    BASALT -a 1_${sample_id}_megahit.fa,2_${sample_id}_spades.fa -s ${sample_id}_1.fq.gz,${sample_id}_2.fq.gz -t 32 -m 64
    
    tar cvfz ${sample_id}.FinalBinSet.tar.gz Final_bestbinset  
    
    """
}


workflow {
    preppedReads = PREP(reads)

    humannProfile = PROFILE_HUMANN(preppedReads)
    metaphlanProfile = PROFILE_METAPHLAN(preppedReads)
    PROFILE_MERGE_HUMANN(humannProfile.collect())
    PROFILE_MERGE_METAPHLAN(metaphlanProfile.collect())

    asmMegahit = ASM_MEGAHIT(preppedReads)
    asmSpades = ASM_SPADES(preppedReads)

    BASALT(preppedReads, asmMegahit, asmSpades)
}
