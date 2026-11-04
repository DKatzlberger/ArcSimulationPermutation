// Pipeline paramters
params {
    data: Path = "data/JB22_Spleen/rna_integrated.rds"

    // Knockout
    knockout = [ "Cysltr2", "Brd7" ] // "Cmklr1", "Cd74", "Arg1", "Shc1"

    // ID of permutation (1 is not shuffled "real")
    seed = [ 1, 42, 34 ]

    // limma
    limma_formula = [ "~genotype" , "~genotype + orig.ident" ]
    limma_ref     = [ "NTC", "NTC" ]

    // edger
    edger_formula = [ "~genotype" , "~genotype + orig.ident" ]
    edger_ref     = [ "NTC", "NTC" ]

    // deseq2
    deseq2_formula = [ "~genotype" , "~genotype + orig.ident" ]
    deseq2_ref     = [ "NTC", "NTC" ]

    // dream
    dream_formula = [ "~genotype" , "~genotype + orig.ident" , "~genotype + (1|orig.ident)" , "~genotype + orig.ident + (1|guide)", "~genotype + (1|orig.ident) + (1|guide)" ]
    dream_ref     = [ "NTC" , "NTC" , "NTC" , "NTC", "NTC" ]

    // nebula
    nebula_formula = [ "~genotype" , "~genotype + orig.ident" , "~genotype + (1|orig.ident)" , "~genotype + orig.ident + (1|guide)" ]
    nebula_ref     = [ "NTC" , "NTC" , "NTC" , "NTC" ]
    
    // sceptre
    sceptre_formula = [ "~genotype" , "~genotype + orig.ident" ]
    sceptre_ref     = [ "NTC" , "NTC" ]
    
    // wilcox
    wilcox_formula = [ "~genotype" ]
    wilcox_ref     = [ "NTC" ]

    // scmageck
    scmageck_formula = [ "~genotype" ]
    scmageck_ref     = [ "NTC" ]
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

process pseudobulk {

    input:
    path data
    val knockout
    val seed

    output:
    tuple path("pseudobulk_knockout_${knockout}_seed_${seed}.h5ad"),
          val(knockout),
          val(seed)

    script:
    """
    Rscript ${workflow.launchDir}/R/pseudobulk.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}"
    """
}

process limma {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("limma_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/limma.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process edger {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("edger_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/edger.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process deseq2 {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("deseq2_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/deseq2.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process dream {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("dream_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/dream.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process nebula {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("nebula_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/nebula.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process sceptre {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("sceptre_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/sceptre.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process wilcox {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("wilcox_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/wilcox.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
    """
}

process scmageck {

    input:
    path data
    val knockout
    val seed
    val formula
    val ref
    
    output:
    path("scmageck_knockout_${knockout}_seed_${seed}_${formula}_dge.rds")
    
    script:
    """
    Rscript ${workflow.launchDir}/R/scmageck.R \
        "${workflow.launchDir}" \
        "${data}" \
        "${knockout}" \
        "${seed}" \
        "${formula}" \
        "${ref}"
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

    // R: Pseudobulk ----------------------------------------------------------
    pseudobulk_ch = shuffle.out

    pseudobulk(
        pseudobulk_ch.map { table -> table[0] }, // data
        pseudobulk_ch.map { table -> table[1] }, // knockout
        pseudobulk_ch.map { table -> table[2] }  // seed
    )

    // DGE pseudobulk =========================================================
    // R: limma ---------------------------------------------------------------
    limma_formula = channel.fromList(params.limma_formula)
    limma_ref     = channel.fromList(params.limma_ref)
    limma_pairs   = limma_formula.merge(limma_ref)

    limma_args = pseudobulk.out.combine(limma_pairs) // [data, knockout, seed, formula, ref]

    limma(
        limma_args.map { table -> table[0] }, // data
        limma_args.map { table -> table[1] }, // knockout
        limma_args.map { table -> table[2] }, // seed
        limma_args.map { table -> table[3] }, // formula
        limma_args.map { table -> table[4] }  // ref
    )

    // R: edger ---------------------------------------------------------------
    edger_formula = channel.fromList(params.edger_formula)
    edger_ref     = channel.fromList(params.edger_ref)
    edger_pairs   = edger_formula.merge(edger_ref)

    edger_args = pseudobulk.out.combine(edger_pairs) // [data, knockout, seed, formula, ref]

    edger(
        edger_args.map { table -> table[0] }, // data
        edger_args.map { table -> table[1] }, // knockout
        edger_args.map { table -> table[2] }, // seed
        edger_args.map { table -> table[3] }, // formula
        edger_args.map { table -> table[4] }  // ref
    )

    // R: deseq2 ---------------------------------------------------------------
    deseq2_formula = channel.fromList(params.deseq2_formula)
    deseq2_ref     = channel.fromList(params.deseq2_ref)
    deseq2_pairs   = deseq2_formula.merge(deseq2_ref)

    deseq2_args = pseudobulk.out.combine(deseq2_pairs) // [data, knockout, seed, formula, ref]

    deseq2(
        deseq2_args.map { table -> table[0] }, // data
        deseq2_args.map { table -> table[1] }, // knockout
        deseq2_args.map { table -> table[2] }, // seed
        deseq2_args.map { table -> table[3] }, // formula
        deseq2_args.map { table -> table[4] }  // ref
    )

    // R: dream ---------------------------------------------------------------
    dream_formula = channel.fromList(params.dream_formula)
    dream_ref     = channel.fromList(params.dream_ref)
    dream_pairs   = dream_formula.merge(dream_ref)

    dream_args = pseudobulk.out.combine(dream_pairs) // [data, knockout, seed, formula, ref]

    dream(
        dream_args.map { table -> table[0] }, // data
        dream_args.map { table -> table[1] }, // knockout
        dream_args.map { table -> table[2] }, // seed
        dream_args.map { table -> table[3] }, // formula
        dream_args.map { table -> table[4] }  // ref
    )

    // DGE single cell ========================================================
    // R: nebula --------------------------------------------------------------
    nebula_formula = channel.fromList(params.nebula_formula)
    nebula_ref     = channel.fromList(params.nebula_ref)
    nebula_pairs   = nebula_formula.merge(nebula_ref)

    nebula_args = shuffle.out.combine(nebula_pairs) // [data, knockout, seed, formula, ref]

    nebula(
        nebula_args.map { table -> table[0] }, // data
        nebula_args.map { table -> table[1] }, // knockout
        nebula_args.map { table -> table[2] }, // seed
        nebula_args.map { table -> table[3] }, // formula
        nebula_args.map { table -> table[4] }  // ref
    )

    // R: sceptre -------------------------------------------------------------
    sceptre_formula = channel.fromList(params.sceptre_formula)
    sceptre_ref     = channel.fromList(params.sceptre_ref)
    sceptre_pairs   = sceptre_formula.merge(sceptre_ref)

    sceptre_args = shuffle.out.combine(sceptre_pairs) // [data, knockout, seed, formula, ref]

    sceptre(
        sceptre_args.map { table -> table[0] }, // data
        sceptre_args.map { table -> table[1] }, // knockout
        sceptre_args.map { table -> table[2] }, // seed
        sceptre_args.map { table -> table[3] }, // formula
        sceptre_args.map { table -> table[4] }  // ref
    )

    // R: wilcox --------------------------------------------------------------
    wilcox_formula = channel.fromList(params.wilcox_formula)
    wilcox_ref     = channel.fromList(params.wilcox_ref)
    wilcox_pairs   = wilcox_formula.merge(wilcox_ref)

    wilcox_args = shuffle.out.combine(wilcox_pairs) // [data, knockout, seed, formula, ref]

    wilcox(
        wilcox_args.map { table -> table[0] }, // data
        wilcox_args.map { table -> table[1] }, // knockout
        wilcox_args.map { table -> table[2] }, // seed
        wilcox_args.map { table -> table[3] }, // formula
        wilcox_args.map { table -> table[4] }  // ref
    )

    // R: scmageck ------------------------------------------------------------
    scmageck_formula = channel.fromList(params.scmageck_formula)
    scmageck_ref     = channel.fromList(params.scmageck_ref)
    scmageck_pairs   = scmageck_formula.merge(scmageck_ref)

    scmageck_args = shuffle.out.combine(scmageck_pairs) // [data, knockout, seed, formula, ref]

    scmageck(
        scmageck_args.map { table -> table[0] }, // data
        scmageck_args.map { table -> table[1] }, // knockout
        scmageck_args.map { table -> table[2] }, // seed
        scmageck_args.map { table -> table[3] }, // formula
        scmageck_args.map { table -> table[4] }  // ref
    )

    publish:
    res_prepare    = prepare.out
    res_shuffle    = shuffle.out
    res_pseudobulk = pseudobulk.out

    res_limma    = limma.out
    res_edger    = edger.out
    res_deseq2   = deseq2.out
    res_dream    = dream.out
    res_nebula   = nebula.out
    res_sceptre  = sceptre.out
    res_wilcox   = wilcox.out
    res_scmageck = scmageck.out
}

// Publish ====================================================================
output {
    res_prepare {
        path "nextflow/ArcSimPerm_JB22_Bcells/prepare"
        mode "copy"
    }
    res_shuffle {
        path "nextflow/ArcSimPerm_JB22_Bcells/shuffle"
        mode "copy"
    }
    res_pseudobulk {
        path "nextflow/ArcSimPerm_JB22_Bcells/pseudobulk"
        mode "copy"
    }

    // DGE results
    res_limma {
        path "nextflow/ArcSimPerm_JB22_Bcells/limma"
        mode "copy"
    }
    res_edger {
        path "nextflow/ArcSimPerm_JB22_Bcells/edger"
        mode "copy"
    }
    res_deseq2 {
        path "nextflow/ArcSimPerm_JB22_Bcells/deseq2"
        mode "copy"
    }
    res_dream {
        path "nextflow/ArcSimPerm_JB22_Bcells/dream"
        mode "copy"
    }
    res_nebula {
        path "nextflow/ArcSimPerm_JB22_Bcells/nebula"
        mode "copy"
    }
    res_sceptre {
        path "nextflow/ArcSimPerm_JB22_Bcells/sceptre"
        mode "copy"
    }
    res_wilcox {
        path "nextflow/ArcSimPerm_JB22_Bcells/wilcox"
        mode "copy"
    }
    res_scmageck {
        path "nextflow/ArcSimPerm_JB22_Bcells/scmageck"
        mode "copy"
    }
}