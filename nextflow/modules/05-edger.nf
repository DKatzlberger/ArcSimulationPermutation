/*
 * Uses the R library 'edgeR' to do DGE analysis.
 */

params.script_id = "05"

 process edger {

    input:
    path data
    val formula
    val ref
    
    output:
    path("${params.script_id}-${data.baseName}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/05-edger.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${formula}" \
        "${ref}"
    """
}