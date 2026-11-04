/*
 * Uses the R library 'scMAGeCK' to do DGE analysis.
 */

params.script_id = "11"

process scmageck {

    input:
    path data
    val formula
    val ref
    
    output:
    path("${params.script_id}-${data.baseName}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/11-scmageck.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${formula}" \
        "${ref}"
    """
}