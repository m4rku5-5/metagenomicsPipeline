#!/bin/nextflow

cleanReads = Channel.fromFilePairs(params.cleanReads, checkIfExists: true)
//cleanReads.view()
basaltBins = Channel.fromPath(params.basaltBins, checkIfExists: true)
//basaltBins.view()
sampleGroups = Channel.fromPath(params.metaData, checkIfExists: true)
                  .splitCsv()
//sampleGroups.view()

process PREP_BINS {
    cpus 2
    time '15min'
    memory '8 GB'
    
    input:
         path(basaltBin)

    output:
         path("*.ren.tar.gz")
    
    script:
    """
    Rscript 02a.prepBins.R ${basaltBin}
    
    """
}

process ANNOT_PROKKA {
    cpus 24
    time '3h'
    memory '32 GB'
    scratch true
    stageOutMode 'move'
	
    input:
         path(renBin)

    output:
         path("*.prokka.tar.gz")
         
    publishDir "PROKKA", mode: 'copy'


    script:
    """
    conda activate ANNOT
        
    sampleID="\$(echo ${renBin} | grep -o 'P[^.]*')"
    
    mkdir PROKKA_\${sampleID}
    cd PROKKA_\${sampleID}
    tar xvfz ../${renBin}
        
    mkdir fixed
    doawk(){ awk '/^>/ { printf("%s_%s\\n",\$0,i++);next;} { print \$0;}' \$1 > fixed/\$2.fa; }
    export -f doawk
    ls \${sampleID}/*.fa | parallel -j4 doawk {} {/.}
    
    ls fixed/*.fa | parallel -j2 prokka {} --outdir out/\${sampleID}-{/.} --locustag \${sampleID}-{/.} --prefix \${sampleID}-{/.} --force
    
    mkdir need
    mv out/*/*.gff need
    mv out/*/*.gbk need
    mv out/*/*.faa need
    mv out/*/*.ffn need
    
    tar czvf ../\${sampleID}.prokka.tar.gz need
    
    """
}

process CLUSTER {
    cpus 32
    time '10h'
    memory '196 GB'
    
    input:
         path("*.prokka.tar.gz")

    output:
         path("proteinCatalogue_nr*")
         
    publishDir "PROTCAT", mode: 'copy'
    
    script:
    """
    conda activate ANNOT    
    
    ls *.prokka.tar.gz | parallel -j1 'tar -axOf {} --wildcards "need/*.faa" ' > proteinCatalogue_r.faa
    mmseqs easy-cluster -c 1 --min-seq-id 0.97 --threads 32 proteinCatalogue_r.faa proteinCatalogue_nr tmp
    rm -r tmp
    diamond makedb -d proteinCatalogue_nr --threads 32 --in proteinCatalogue_nr_rep_seq.fasta
    """
}

process ANNOT_EGGNOG {
    cpus 32
    time '16h'
    memory '128 GB'
    
    input:
         path(protcat_nr)

    output:
         path("emapper*")
         
    publishDir "PROTCAT", mode: 'copy'
    
    script:
    """
    conda activate ANNOT
    
    emapper.py --data_dir DB/EMData -i proteinCatalogue_nr_rep_seq.fasta -o emapper --cpu 32
    """
}

process GROUP_DREP {
    cpus 24
    time '6h'
    memory '64 GB'
	scratch true
    stageOutMode 'move'
    
    input:
         tuple path(renBins), val(donor) 

    output:
         path("derepBins-${donor}.tar.gz")
         
    publishDir "BINCAT", mode: 'copy'
    
    script:
    """
    conda activate BINPROC
    
    ls *.ren.tar.gz | xargs -i tar xvfz {}
    
    echo "genome,contamination,completeness" > genomeInfo.csv
    awk FNR!=1 P*/bins_checkm.csv >> genomeInfo.csv
    
    mkdir bins
    mv P*/*.fa bins
    
    mkdir ${donor}
    mv bins ${donor}
    mv genomeInfo.csv ${donor}
    cd ${donor}
    
    dRep dereplicate dereped-${donor} -g bins/*.fa --genomeInfo genomeInfo.csv -p 24
    
    tar czvf ../derepBins-${donor}.tar.gz dereped-${donor}
        
    """
}

process COMPLETE_DREP {
    cpus 24
    time '6h'
    memory '64 GB'
	scratch true
    stageOutMode 'move'
    
    input:
         path(renBins)

    output:
         path("derepBins-complete.tar.gz")
         
    publishDir "BINCAT", mode: 'copy'
    
    script:
    """
    conda activate BINPROC

    mkdir COMPLETEDEREP
    mv *.ren.tar.gz COMPLETEDEREP
    cd COMPLETEDEREP
    
    ls *.ren.tar.gz | xargs -i tar xvfz {}
    
    echo "genome,contamination,completeness" > genomeInfo.csv
    awk FNR!=1 P*/bins_checkm.csv >> genomeInfo.csv
    
    mkdir bins
    mv P*/*.fa bins 
    
    dRep dereplicate derepBins-COMPLETE -g bins/*.fa --genomeInfo genomeInfo.csv -p 24
    
    tar czvf ../derepBins-complete.tar.gz derepBins-COMPLETE
    
    """
}


process MAP_MAKEINDEX {
    cpus 24
    time '6h'
    memory '128 GB'
    
    input:
         path(derepBins_complete)

    output:
         path("BINCAT.index*")
         
    publishDir "MAP", mode: 'copy'
    
    script:
    """
    conda activate BINPROC    
    
    tar -axOf ${derepBins_complete} --wildcards "derepBins-COMPLETE/dereplicated_genomes/*.fa" | awk '/^>/ { printf("%s_%s\\n",\$0,i++);next;} { print \$0;}' > BINCAT.fna
    
    bowtie2-build BINCAT.fna BINCAT.index --threads 24
    
    rm BINCAT.fna
    
    """
}

process MAP_MAP {
    cpus 24
    time '3h'
    memory '128 GB'
    
    input:
         tuple val(sample_id), path(reads)
         path(binIndex)

    output:
         path("${sample_id}.s.bam")
         
    publishDir "MAP", mode: 'symlink'
    
    script:
    """
    conda activate BINPROC
    
    idx=\$(echo ${binIndex[0]} | sed 's#.1.bt2##g')
    
    bowtie2 -p 24 -x \${idx} -1 ${reads[0]} -2 ${reads[1]} | samtools view -@ 4 -b -F 0x900 -q 1 - > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}.s.bam

    """
}

process MAP_ABUND {
    cpus 4
    time '30min'
    memory '12 GB'
    
    input:
         path(mappedBam)

    output:
         path("*.idxstat")
         
    publishDir "BINABUND", mode: 'copy'
    
    script:
    """
    conda activate BINPROC
    
    sampleID=\$(basename -s .bam ${mappedBam})
    
    samtools idxstats ${mappedBam} > \${sampleID}.idxstat

    """
}

process TAX {
    cpus 24
    time '3h'
    memory '250 GB'
    
    input:
         path(derepBins_complete)

    output:
         path("binTax/*")
         
    publishDir "TAX", mode: 'copy'
    
    script:
    """
    conda activate GTDBTK

    tar xvfz ${derepBins_complete} --wildcards "derepBins-COMPLETE/dereplicated_genomes/*.fa"

    gtdbtk classify_wf --genome_dir derepBins-COMPLETE/dereplicated_genomes/ --extension .fa --out_dir binTax --cpus 24 --pplacer_cpus 12
         
    """
}



workflow {
   preppedBins = PREP_BINS(basaltBins)
   prokka = ANNOT_PROKKA(preppedBins)
   protcat = CLUSTER(prokka.collect())
   ANNOT_EGGNOG(protcat)
   mergedGroups = preppedBins.map{x -> [x.name.replaceFirst(~/\.[^\.]+$/, '').split("\\.")[0], x]}
    .join(sampleGroups)
    .groupTuple(by: 2)
    .map{x -> [x[1], x[2]]}

   GROUP_DREP(mergedGroups)
   completeDrep = COMPLETE_DREP(preppedBins.collect())
   TAX(completeDrep)
   
   binIndex = MAP_MAKEINDEX(completeDrep)
   mappedBam = MAP_MAP(cleanReads, binIndex)
   MAP_ABUND(mappedBam)
}
