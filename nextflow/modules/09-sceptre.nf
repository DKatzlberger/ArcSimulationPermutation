/*
 * Uses the R library 'sceptre' to do DGE analysis.
 */

params.script_id = "09"

process sceptre {

    input:
    path data
    val formula
    val ref
    
    output:
    path("${params.script_id}-${data.baseName}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/09-sceptre.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${formula}" \
        "${ref}"
    """
}