#!/usr/bin/env nextflow
/****************************************
 * Metagenomics profiling workflow 
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 ****************************************/
version = '0.9b1'

// Create file objects
kaiju_db = file(params.kaiju_db)
kaiju_nodes = file(params.kaiju_nodes)
kaiju_names = file(params.kaiju_names)

// Channels with paired input reads
Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_kaiju;
           input_reads_bbmap}

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

process summarize_taxonomic_profile {
    tag {merged_taxonomic_profiles.baseName}
    publishDir "${params.outdir}/taxonomic_profiles", mode: 'copy'

    input:
    file taxonomic_profiles from kaiju_output
    file nodes from kaiju_nodes
    file names from kaiju_names

    output:
    file "${taxonomic_profiles.baseName}.krona"
    file "${taxonomic_profiles.baseName}.krona.html"
    file "${taxonomic_profiles.baseName}.summary.{family,genus,species}"

    """
    kaiju2krona \
        -t $nodes \
        -n $names \
        -i $taxonomic_profiles \
        -o ${taxonomic_profiles.baseName}.krona
    ktImportText \
        -o ${taxonomic_profiles.baseName}.krona.html \
        ${taxonomic_profiles.baseName}.krona
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

