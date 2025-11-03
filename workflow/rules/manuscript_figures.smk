rule plot_figure_1d:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_gamma=f"{OUTPUTS_DIR}/gi_scores/screen2023/construct_scores/all_gis_Gamma.OI.Avg.tsv"
    output:
        output_figure_1d=f"{OUTPUTS_DIR}/manuscript_figures/figure_1d.svg"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_1d.log"
    script:
        "../scripts/manuscript_figures/plot_figure_1d.R"

rule plot_figure_1e:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_tau=f"{OUTPUTS_DIR}/gi_scores/screen2023/construct_scores/all_gis_Tau.OI.Avg.tsv"
    output:
        output_figure_1e=f"{OUTPUTS_DIR}/manuscript_figures/figure_1e.svg"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_1e.log"
    script:
        "../scripts/manuscript_figures/plot_figure_1e.R"

rule plot_figure_2a:
    input:
        # updated paths: differential outputs are now placed with Nu.* naming under construct/gene/discriminant dirs
        input_nu_r1=f"{OUTPUTS_DIR}/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R1.tsv",
        input_nu_r2=f"{OUTPUTS_DIR}/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R2.tsv",
        input_nu_gene_level=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv"
    output:
        output_figure_2a=f"{OUTPUTS_DIR}/manuscript_figures/figure_2a.png",
        output_figure_2a_labels=f"{OUTPUTS_DIR}/manuscript_figures/figure_2a_labels.png"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_2a.log"
    script:
        "../scripts/manuscript_figures/plot_figure_2a.R"

rule plot_figure_2b:
    input:
        input_nu_r1=f"{OUTPUTS_DIR}/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R1.tsv",
        input_nu_r2=f"{OUTPUTS_DIR}/gi_scores/screen2023/gene_combination_scores/gene_combination_scores_Nu.OI.R2.tsv",
        input_nu_gene_level=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_idmap=(
            "data/annotations/screen2023_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "counts_metadata" 
            else f"{OUTPUTS_DIR}/annotations/screen2023_gene_combination_id_map.tsv"
        )
    output:
        output_figure_2b=f"{OUTPUTS_DIR}/manuscript_figures/figure_2b.png",
        output_figure_2b_labels=f"{OUTPUTS_DIR}/manuscript_figures/figure_2b_labels.png"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_2b.log"
    script:
        "../scripts/manuscript_figures/plot_figure_2b.R"

rule plot_figure_2d:
    input:
        input_nu=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_idmap=(
            "data/annotations/screen2023_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "counts_metadata" 
            else f"{OUTPUTS_DIR}/annotations/screen2023_gene_combination_id_map.tsv"
        )
    output:
        output_figure_2d_svg=f"{OUTPUTS_DIR}/manuscript_figures/figure_2d.svg"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_2d.log"
    script:
        "../scripts/manuscript_figures/plot_figure_2d.R"

rule plot_figure_3a:
    input:
        input_gamma=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_idmap=(
            "data/annotations/screen2023_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "counts_metadata" 
            else f"{OUTPUTS_DIR}/annotations/screen2023_gene_combination_id_map.tsv"
        )
    output:
        output_figure_3a=f"{OUTPUTS_DIR}/manuscript_figures/figure_3a.png"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_3a.log"
    script:
        "../scripts/manuscript_figures/plot_figure_3a.R"

rule plot_figure_3b_negative:
    input:
        input_gamma=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_idmap=(
            "data/annotations/screen2023_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "counts_metadata" 
            else f"{OUTPUTS_DIR}/annotations/screen2023_gene_combination_id_map.tsv"
        )
    output:
        output_figure_3b=f"{OUTPUTS_DIR}/manuscript_figures/figure_3b_negative.png"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_3b_negative.log"
    script:
        "../scripts/manuscript_figures/plot_figure_3b_negative.R"


rule plot_figure_3b_positive:
    input:
        input_gamma=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Gamma.OI.Avg.tsv",
        input_tau=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Tau.OI.Avg.tsv",
        input_nu=f"{OUTPUTS_DIR}/gi_scores/screen2023/discriminant_scores/discriminant_hits_Nu.OI.Avg.tsv",
        input_idmap=(
            "data/annotations/screen2023_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "counts_metadata" 
            else f"{OUTPUTS_DIR}/annotations/screen2023_gene_combination_id_map.tsv"
        )
    output:
        output_figure_3b=f"{OUTPUTS_DIR}/manuscript_figures/figure_3b_positive.png"
    log:
        f"{LOGS_DIR}/manuscript_figures/plot_figure_3b_positive.log"
    script:
        "../scripts/manuscript_figures/plot_figure_3b_positive.R"