// Modules
include { prepare }    from "../modules/01-prepare.nf"
include { shuffle }    from "../modules/02-shuffle.nf"
include { pseudobulk } from "../modules/03-pseudobulk.nf"
include { limma }      from "../modules/04-limma.nf"
include { edger }      from "../modules/05-edger.nf"
include { deseq2 }     from "../modules/06-deseq2.nf"
include { dream }      from "../modules/07-dream.nf"
include { nebula }     from "../modules/08-nebula.nf"
include { sceptre }    from "../modules/09-sceptre.nf"
include { wilcox }     from "../modules/10-wilcox.nf"
include { scmageck }   from "../modules/11-scmageck.nf"
include { visualize }  from "../modules/12-visualize.nf"

// Paramters
params {
    data: Path = "data/GSE213511/rna_integrated.rds"

    // Prepare
    prepare_celltype = [ "Ery", "Imm. B-cell", "Gran. P", "GMP", "Gran.", "HSC" ] // "Ery", "Imm. B-cell", "Gran. P", "GMP", "Gran.", "HSC"
    prepare_knockout = [ "Cebpa", "Kmt2d", "Brd9", "Rcor1" ] // "Cebpa", "Kmt2d", "Brd9", "Rcor1"

    // Shuffle (1 is not shuffled "real")
    shuffle_seed = (1..100)

    // limma
    limma_formula = [ "~genotype" , "~genotype + sample" ]
    limma_ref     = [ "NTC", "NTC" ]

    // edger
    edger_formula = [ "~genotype" , "~genotype + sample" ]
    edger_ref     = [ "NTC", "NTC" ]

    // deseq2
    deseq2_formula = [ "~genotype" , "~genotype + sample" ]
    deseq2_ref     = [ "NTC", "NTC" ]

    // dream
    dream_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)", "~genotype + (1|sample) + (1|guide)" ]
    dream_ref     = [ "NTC" , "NTC" , "NTC" , "NTC", "NTC" ]

    // nebula
    nebula_formula = [ "~genotype" , "~genotype + sample" , "~genotype + (1|sample)" , "~genotype + sample + (1|guide)" ]
    nebula_ref     = [ "NTC" , "NTC" , "NTC" , "NTC" ]
    
    // sceptre
    sceptre_formula = [ "~genotype" , "~genotype + sample" ]
    sceptre_ref     = [ "NTC" , "NTC" ]
    
    // wilcox
    wilcox_formula = [ "~genotype" ]
    wilcox_ref     = [ "NTC" ]

    // scmageck
    scmageck_formula = [ "~genotype" ]
    scmageck_ref     = [ "NTC" ]
}


// Workflow ===================================================================
workflow {
    main:

    // R: prepare -------------------------------------------------------------
    prepare_args = channel.value(params.data)
        .combine(channel.fromList(params.prepare_celltype))
        .combine(channel.fromList(params.prepare_knockout))

    prepare(
        prepare_args.map { arg -> arg[0] }, // data
        prepare_args.map { arg -> arg[1] }, // celltype
        prepare_args.map { arg -> arg[2] }  // knockout
    )

    // R: shuffle -------------------------------------------------------------
    shuffle_args = prepare.out.combine(channel.fromList(params.shuffle_seed))

    shuffle(
        shuffle_args.map { arg -> arg[0] }, // data
        shuffle_args.map { arg -> arg[1] }  // seed
    )

    // R: pseudobulk ----------------------------------------------------------
    pseudobulk_args = shuffle.out

    pseudobulk(
        pseudobulk_args
    )

   // DGE pseudobulk =========================================================
    // R: limma ---------------------------------------------------------------
    limma_formula = channel.fromList(params.limma_formula)
    limma_ref     = channel.fromList(params.limma_ref)
    limma_pairs   = limma_formula.merge(limma_ref)

    limma_args = pseudobulk.out.combine(limma_pairs)

    limma(
        limma_args.map { arg -> arg[0] }, // data
        limma_args.map { arg -> arg[1] }, // formula
        limma_args.map { arg -> arg[2] }  // ref
    )

    // R: edger ---------------------------------------------------------------
    edger_formula = channel.fromList(params.edger_formula)
    edger_ref     = channel.fromList(params.edger_ref)
    edger_pairs   = edger_formula.merge(edger_ref)

    edger_args = pseudobulk.out.combine(edger_pairs)

    edger(
        edger_args.map { arg -> arg[0] }, // data
        edger_args.map { arg -> arg[1] }, // formula
        edger_args.map { arg -> arg[2] }  // ref
    )

    // R: deseq2 --------------------------------------------------------------
    deseq2_formula = channel.fromList(params.deseq2_formula)
    deseq2_ref     = channel.fromList(params.deseq2_ref)
    deseq2_pairs   = deseq2_formula.merge(deseq2_ref)

    deseq2_args = pseudobulk.out.combine(deseq2_pairs)

    deseq2(
        deseq2_args.map { arg -> arg[0] }, // data
        deseq2_args.map { arg -> arg[1] }, // formula
        deseq2_args.map { arg -> arg[2] }  // ref
    )

    // R: dream ---------------------------------------------------------------
    dream_formula = channel.fromList(params.dream_formula)
    dream_ref     = channel.fromList(params.dream_ref)
    dream_pairs   = dream_formula.merge(dream_ref)

    dream_args = pseudobulk.out.combine(dream_pairs)

    dream(
        dream_args.map { arg -> arg[0] }, // data
        dream_args.map { arg -> arg[1] }, // formula
        dream_args.map { arg -> arg[2] }  // ref
    )

    // DGE single cell ========================================================
    // R: nebula --------------------------------------------------------------
    nebula_formula = channel.fromList(params.nebula_formula)
    nebula_ref     = channel.fromList(params.nebula_ref)
    nebula_pairs   = nebula_formula.merge(nebula_ref)

    nebula_args = shuffle.out.combine(nebula_pairs)

    nebula(
        nebula_args.map { arg -> arg[0] }, // data
        nebula_args.map { arg -> arg[1] }, // formula
        nebula_args.map { arg -> arg[2] }  // ref
    )

    // R: sceptre -------------------------------------------------------------
    sceptre_formula = channel.fromList(params.sceptre_formula)
    sceptre_ref     = channel.fromList(params.sceptre_ref)
    sceptre_pairs   = sceptre_formula.merge(sceptre_ref)

    sceptre_args = shuffle.out.combine(sceptre_pairs)

    sceptre(
        sceptre_args.map { arg -> arg[0] }, // data
        sceptre_args.map { arg -> arg[1] }, // formula
        sceptre_args.map { arg -> arg[2] }  // ref
    )

    // R: wilcox --------------------------------------------------------------
    wilcox_formula = channel.fromList(params.wilcox_formula)
    wilcox_ref     = channel.fromList(params.wilcox_ref)
    wilcox_pairs   = wilcox_formula.merge(wilcox_ref)

    wilcox_args = shuffle.out.combine(wilcox_pairs)

    wilcox(
        wilcox_args.map { arg -> arg[0] }, // data
        wilcox_args.map { arg -> arg[1] }, // formula
        wilcox_args.map { arg -> arg[2] }  // ref
    )

    // R: scmageck ------------------------------------------------------------
    scmageck_formula = channel.fromList(params.scmageck_formula)
    scmageck_ref     = channel.fromList(params.scmageck_ref)
    scmageck_pairs   = scmageck_formula.merge(scmageck_ref)

    scmageck_args = shuffle.out.combine(scmageck_pairs)

    scmageck(
        scmageck_args.map { arg -> arg[0] }, // data
        scmageck_args.map { arg -> arg[1] }, // formula
        scmageck_args.map { arg -> arg[2] }  // ref
    )

    // Visualize ==============================================================
    // R: 15-visualize --------------------------------------------------------
    dge_results = limma.out.collect()
        .mix(edger.out)
        .mix(deseq2.out)
        .mix(dream.out)
        .mix(nebula.out)
        .mix(sceptre.out)
        .mix(wilcox.out)
        .mix(scmageck.out)
        .collect()

    visualize(
        dge_results
    )

    publish:
    // Shuffle res
    res_prepare    = prepare.out
    res_shuffle    = shuffle.out
    res_pseudobulk = pseudobulk.out

    // DGE res
    res_limma    = limma.out
    res_edger    = edger.out
    res_deseq2   = deseq2.out
    res_dream    = dream.out
    res_nebula   = nebula.out
    res_sceptre  = sceptre.out
    res_wilcox   = wilcox.out
    res_scmageck = scmageck.out

    // Visualize
    res_visualize = visualize.out.panels
}

// Publish ====================================================================
output {
    // Shuffle res
    res_prepare {
        path "01-prepare"
        mode "copy"
    }
    res_shuffle {
        path "02-shuffle"
        mode "copy"
    }
    res_pseudobulk {
        path "03-pseudobulk"
        mode "copy"
    }

    // DGE res
    res_limma {
        path "04-limma"
        mode "copy"
    }
    res_edger {
        path "05-edger"
        mode "copy"
    }
    res_deseq2 {
        path "06-deseq2"
        mode "copy"
    }
    res_dream {
        path "07-dream"
        mode "copy"
    }
    res_nebula {
        path "08-nebula"
        mode "copy"
    }
    res_sceptre {
        path "09-sceptre"
        mode "copy"
    }
    res_wilcox {
        path "10-wilcox"
        mode "copy"
    }
    res_scmageck {
        path "11-scmageck"
        mode "copy"
    }

    // Visualize
    res_visualize {
        path "12-visualize"
        mode "copy"
    }
}