#!/usr/bin/env nextflow
/****************************************
 * Prototype CTMR metagenomics workflow
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 ****************************************/
version = '0.1a'

// Create file objects
kaiju_db = file(params.kaiju_db)
kaiju_nodes = file(params.kaiju_nodes)
kaiju_names = file(params.kaiju_names)
bracken_kmer_distribution = file(params.bracken_kmer_distribution)

// Channels with paired input reads
Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_kaiju;
           input_reads_kraken;
           input_reads_bbmap;
           input_reads_paladin}

/****************************************
 *    TAXONOMIC COMPOSITION ESTIMATION
 ****************************************/
process kaiju {
    tag {pair_id}
    publishDir "${params.outdir}/kaiju", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_kaiju
    file database from kaiju_db
    file nodes from kaiju_nodes

    output:
    file "${pair_id}.kaiju" into kaiju_output

    """
    if [[ "${reads[0]}" == *.gz ]]
    then
        kaiju \
            -z ${task.cpus} \
            -t $nodes \
            -f $database \
            -i <(gunzip -c ${reads[0]}) \
            -j <(gunzip -c ${reads[1]}) \
            -o ${pair_id}.kaiju
    else
        kaiju \
            -z ${task.cpus} \
            -t $nodes \
            -f $database \
            -i ${reads[0]} \
            -j ${reads[1]} \
            -o ${pair_id}.kaiju
    fi
    """
} 


process kraken {
    tag {pair_id}
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_kraken

    output:
    file "${pair_id}.kraken" into kraken_output
    file "${pair_id}.kraken.report" into kraken_report

    """
    if [[ "${reads[0]}" == *.gz ]]
    then
        kraken \
            --preload \
            --db ${params.kraken_db} \
            --threads ${task.cpus} \
            --paired \
            --gzip-compressed \
            --fastq-input \
            ${reads[0]} \
            ${reads[1]} \
            --output ${pair_id}.kraken
    else
        kraken \
            --preload \
            --db ${params.kraken_db} \
            --threads ${task.cpus} \
            --paired \
            ${reads[0]} \
            ${reads[1]} \
            --output ${pair_id}.kraken
    fi
    kraken-report \
        --db ${params.kraken_db} \
        --show-zeros \
        ${pair_id}.kraken \
        > ${pair_id}.kraken.report
    """
}


process bracken {
    tag {report}
    publishDir "${params.outdir}/bracken", mode: 'copy'
    
    input:
    file report from kraken_report
    file "kmer_distribution.txt" from bracken_kmer_distribution
    
    output:
    file "${report.baseName}.bracken"
    
    """
    bracken_estimate_abundance.py \
        -i ${report} \
        -k kmer_distribution.txt \
        -l ${params.bracken_classification_level} \
        -t ${params.bracken_threshold} \
        -o ${report.baseName}.bracken
    """
}


process merge_kaiju_kraken {
    tag {kaiju_taxcomp.baseName}
    publishDir "${params.outdir}/merged_taxonomic_classifications", mode: 'copy'

    input:
    file kaiju_taxcomp from kaiju_output
    file kraken_taxcomp from kraken_output

    output:
    file "${kaiju_taxcomp.baseName}.merged_classifications.tab" into merged_kaiju_kraken

    """
    merge_taxonomic_classifications.py \
        --kaiju $kaiju_taxcomp \
        --kraken $kraken_taxcomp \
        --merge-order ${params.taxonomic_merge_order} \
        --dbfile merged_classifications.sqlite3 \
        --output ${kaiju_taxcomp.baseName}.merged_classifications.tab 
    """
}


process summarize_taxonomic_profile {
    tag {merged_taxonomic_profiles.baseName}
    publishDir "${params.outdir}/taxonomic_profiles", mode: 'copy'

    input:
    file merged_taxonomic_profiles from merged_kaiju_kraken
    file nodes from kaiju_nodes
    file names from kaiju_names

    output:
    file "${merged_taxonomic_profiles.baseName}.krona"
    file "${merged_taxonomic_profiles.baseName}.krona.html"
    file "${merged_taxonomic_profiles.baseName}.summary.{family,genus,species}"

    """
    kaiju2krona \
        -t $nodes \
        -n $names \
        -i $merged_taxonomic_profiles \
        -o ${merged_taxonomic_profiles.baseName}.krona
    ktImportText \
        -o ${merged_taxonomic_profiles.baseName}.krona.html \
        ${merged_taxonomic_profiles.baseName}.krona
    for rank in family genus species; do
        kaijuReport \
            -t $nodes \
            -n $names \
            -i $merged_taxonomic_profiles \
            -r \$rank \
            -o ${merged_taxonomic_profiles.baseName}.summary.\$rank
    done
    """
}

process bbmap_igc {
    tag {pair_id}
    publishDir "${params.outdir}/IGC", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_bbmap

    output:
    file "${pair_id}.scafstats.txt.gz" 
    file "${pair_id}.mapping_stats.txt.gz"
    file "${pair_id}.rpkm.txt.gz"
    file "${pair_id}.covstats.txt.gz"
    file "${pair_id}.sam.gz"

    """
    bbmap.sh \
        minid=${params.bbmap_minid} \
        threads=${task.cpus} \
        path=${params.bbmap_igc_dir} \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out=${pair_id}.sam.gz \
        scafstats=${pair_id}.scafstats.txt.gz \
        statsfile=${pair_id}.mapping_stats.txt.gz \
		covstats=${pair_id}.covstats.txt.gz \
		rpkm=${pair_id}.rpkm.txt.gz \
    """
}


process paladin_uniprot {
    tag {pair_id}
    publishDir "${params.outdir}/paladin", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_paladin

    output:
    file "${pair_id}.sam"
    file "${pair_id}_uniprot.tsv"

    """
    zcat ${reads} > all_reads.fastq
    paladin align \
        -t ${task.cpus} \
        -o ${pair_id} \
        ${params.paladin_db} \
        all_reads.fastq
    """
}
