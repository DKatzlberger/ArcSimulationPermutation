params {
    data: Path = "data/JB22_Spleen/rna_integrated.rds"

    // Knockouts of interest
    knockout     = [ "Cysltr2", "Brd7", "Cmklr1", "Cd74", "Arg1", "Shc1" ]

    // ID of permutation
    shuffle_seed = [ 42, 34, 334, 5656, 23432, 3344, 3435, 35345, 4545, 34535 ]

}

// Process ====================================================================
process bpcells {

    input:
    path data
    val  experiment

    output:
   tuple path("bpcells_${experiment}_RNA"), 
         path("bpcells_${experiment}_CRISPR"), 
         path("bpcells_${experiment}.rds")

    script:
    """
    Rscript ${workflow.launchDir}/R/bpcells.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${experiment}"
    """
}

process shuffle {

    input:
    path data
    val seed

    output:
    tuple val(seed),
          path("shuffle_seed_${seed}.rds")


    script:
    """
    Rscript ${workflow.launchDir}/R/shuffle.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${seed}"
    """
}

// Workflow ===================================================================
workflow {
    main:

    // Real data ==============================================================
    bpcells(
        channel.value(params.data),                  // data
        channel.value(file(params.data).parent.name) // experiment
    )

    // Shuffle data ===========================================================

    publish:
    res_bpcells = bpcells.out
}

// Publish ====================================================================
output {
    res_bpcells {
        path "nextflow/bpcells"
        mode "copy"
    }
}