/*
 * Takes a SingleCellExperiment object and prepares it for simulation permutation. It shuffles knockout labels.
 */

params.script_id = "02"

process shuffle {

    input:
    path data
    val seed

    output:
    path("${params.script_id}-${data.baseName}_seed_${seed}.h5ad")

    script:
    """
    Rscript ${workflow.launchDir}/R/02-shuffle.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${seed}"
    """
}