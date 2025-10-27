rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_filtered_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_single_sgRNA_phenotypes.tsv"
    output:
        output_gis="../outputs/gis/{screen}_genetic_interaction_scores.tsv"
    params:
        scores=lambda wildcards: config[f"{wildcards.screen.upper()}_GI_SCORES"]
    script:
        "../scripts/calculate_gi_scores.R"