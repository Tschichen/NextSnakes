FQSAMPLES = null
TRSAMPLES = null
MAMSAMPLES = null

if (PAIRED == 'paired'){
    FR1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/FASTQ/"+element+"_R1.fastq.gz"
    }
    FR2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/FASTQ/"+element+"_R2.fastq.gz"
    }
    FQSAMPLES = FR1+FR2
    FQSAMPLES.sort()

    TR1 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
    }
    TR2 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
    }
    TRSAMPLES = TR1+TR2
    TRSAMPLES.sort()

}else{
    FQSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/FASTQ/"+element+".fastq.gz"
    }
    FQSAMPLES.sort()

    TRSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
    }
    TRSAMPLES.sort()
}

MRSAMPLES = LONGSAMPLES.collect{
    element -> return "${workflow.workDir}/MAPPED/"+element+"_mapped_sorted.sam.gz"
}
MRSAMPLES.sort()

process qc_raw{
    conda 'nextsnakes/envs/qc.yaml'
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

process qc_trimmed{
    conda 'nextsnakes/envs/qc.yaml'
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path read

    output:
    path "*.{zip,html}", emit: trfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

process qc_mapped{
    conda 'nextsnakes/envs/qc.yaml'
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path mapped

    output:
    path "*.{zip,html}", emit: mapfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $read
    """
}


//collecting list of processed file for multiqc, not implemented yet
process collect_qc_raw{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_trimmed{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_map{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

process multiqc{
    conda 'nextsnakes/envs/qc.yaml'
    cpus THREADS
    validExitStatus 0,1
    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$CONDITION/$filename"
        else null
    }

    input:
    path qcs
    path trimmed
    path mapped
    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/QC/FASTQC/$CONDITION/.
    """
}

workflow {
    samples_ch = Channel.from(FQSAMPLES)
    trsamples_ch = Channel.from(TRSAMPLES)
    mapsamples_ch = Channel.from(MAPSAMPLES)

    main:
    qc_raw(samples_ch)
    qc_trimmed(trsamples_ch)
    qc_mapped(mapsamples_ch)

    multiqc_raw(qc_raw.out.fastqc_results, qc_trimmed.out.trfastqc_results, qc_mapped.out.mapfastqc_results)
    emit:
    multiqc.out.multiqc_results
    //collect_qc_raw()
    //multiqc_raw(collect_qc_raw.out.collect_fastqc)
}
