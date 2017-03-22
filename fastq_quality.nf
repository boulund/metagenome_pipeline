#!/usr/bin/env nextflow
/****************************************
 * CTMR metagenomics read QA workflow
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 ****************************************/
version = '0.2a'

params.input_reads = '' // Specify on command line
params.outdir = '.'
params.bbduk_minlen = 75
params.bbduk_qtrim = 'rl'
params.bbduk_trimq = 20
params.bbduk_ktrim = 'r'
params.bbduk_k = 25
params.bbduk_mink = 11
params.bbduk_adapters = '/home/ctmr/src/bbmap/resources/adapters.fa'
params.bbduk_hdist = 1

Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_bbduk}

process bbduk {
    tag {pair_id}
    publishDir "${params.outdir}/qa_reads", mode: 'copy'

    cpus 4
    memory 32.GB
    time 20.m

    input:
    set pair_id, file(reads) from input_reads_bbduk

    output:
    file "${pair_id}_1.fq.gz" 
    file "${pair_id}_2.fq.gz" 
    file "${pair_id}.stats.txt.gz"
    file "${pair_id}.bhist.txt.gz"
    file "${pair_id}.qhist.txt.gz"
    file "${pair_id}.qchist.txt.gz"
    file "${pair_id}.aqhist.txt.gz"
    file "${pair_id}.bqhist.txt.gz"
    file "${pair_id}.lhist.txt.gz"
    file "${pair_id}.gchist.txt.gz"

    """
    bbduk.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        ref=${params.bbduk_ref} \
        out1=${pair_id}_1.fq.gz \
        out2=${pair_id}_2.fq.gz \
        stats=${pair_id}.stats.txt.gz \
        bhist=${pair_id}.bhist.txt.gz \
        qhist=${pair_id}.qhist.txt.gz \
        qchist=${pair_id}.qchist.txt.gz \
        aqhist=${pair_id}.aqhist.txt.gz \
        bqhist=${pair_id}.bqhist.txt.gz \
        lhist=${pair_id}.lhist.txt.gz \
        gchist=${pair_id}.gchist.txt.gz \
        minlen=${params.bbduk_minlen} \
        qtrim=${params.bbduk_qtrim} \
        trimq=${params.bbduk_trimq} \
        ktrim=${params.bbduk_ktrim} \
        k=${params.bbduk_k} \
        mink=${params.bbduk_mink} \
        hdist=${params.bbduk_hdist} \
    """
} 

