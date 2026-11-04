/*
 * Uses the R library 'dream' to do DGE analysis.
 */

params.script_id = "07"

process dream {

    input:
    path data
    val formula
    val ref
    
    output:
    path("${params.script_id}-${data.baseName}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/07-dream.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${formula}" \
        "${ref}"
    """
}