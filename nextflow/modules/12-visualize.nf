params.script_id = "12"

process visualize {

    input:
    path dge_files

    output:
    path "*.svg", emit: svg
    path "*.rds", emit: rds


    script:
    """
    Rscript ${workflow.launchDir}/R/12-visualize.R \
        "${workflow.launchDir}" \
        ${dge_files}
    """
}