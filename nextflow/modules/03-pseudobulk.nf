/*
 * Uses the R library 'scuttle' to aggregate scRNA-seq to pseudobulk by sample and guide information.
 */

params.script_id = "03"

process pseudobulk {
    
    input:
    path data
    
    output:
    path("${params.script_id}-${data.baseName}.h5ad")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/03-pseudobulk.R \
        "${workflow.launchDir}" \
        "${data}"
    """
}