// Pipeline paramters
params {
    data: Path = "data/JB22_Spleen/rna_integrated.rds"

    // Knockouts of interest
    knockout = [ "Cysltr2", "Brd7", "Cmklr1", "Cd74", "Arg1", "Shc1" ]

    // ID of permutation
    seed = [ 1, 42, 34, 334, 5656, 23432, 3344, 3435, 35345, 4545, 34535 ]

}

// Process ====================================================================
process prepare {

    input:
    path data
    val knockout

    output:
    tuple path("prepare_knockout_${knockout}.h5ad"),
          val(knockout)

    script:
    """
    Rscript ${workflow.launchDir}/R/prepare.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}"
    """
}

process shuffle {

    input:
    path data
    val knockout
    val seed

    output:
    tuple path("shuffle_knockout_${knockout}_seed_${seed}.h5ad"),
          val(knockout),
          val(seed)

    script:
    """
    Rscript ${workflow.launchDir}/R/shuffle.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}"
    """
}

// Workflow ===================================================================
workflow {
    main:

    // R: Prepare -------------------------------------------------------------
    prepare_ch = channel.value(params.data).combine(channel.fromList(params.knockout)) // [data, knockout]

    prepare(
        prepare_ch.map { table -> table[0] }, // data
        prepare_ch.map { table -> table[1] }  // knockout
    )

    // R: Shuffle -------------------------------------------------------------
    seeds_ch = channel.fromList(params.seed)

    shuffle_ch = prepare.out.combine(seeds_ch) // [data, knockout, seed]

    shuffle(
        shuffle_ch.map { table -> table[0] }, // data
        shuffle_ch.map { table -> table[1] }, // knockout
        shuffle_ch.map { table -> table[2] }  // seed
    )

    publish:
    res_prepare = prepare.out
    res_shuffle = shuffle.out

}

// Publish ====================================================================
output {
    res_prepare {
        path "nextflow/prepare"
        mode "copy"
    }
    res_shuffle {
        path "nextflow/shuffle"
        mode "copy"
    }
}