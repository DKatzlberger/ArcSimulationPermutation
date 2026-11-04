/*
 * Uses the R library 'limma' to do DGE analysis.
 */

params.script_id = "04"

process limma {

    input:
    path data
    val formula
    val ref
    
    output:
    path("${params.script_id}-${data.baseName}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/04-limma.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${formula}" \
        "${ref}"
    """
}