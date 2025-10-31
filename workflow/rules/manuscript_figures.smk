rule plot_figure_1d:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_gamma="../outputs/gi_scores/screen2023/construct_scores/all_gis_Gamma.OI.Avg.tsv"
    output:
        output_figure_1d="../outputs/manuscript_figures/figure_1d.svg"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_1d.log"
    script:
        "../scripts/manuscript_figures/plot_figure_1d.R"

rule plot_figure_1e:
    input:
        # use the consolidated construct scores file and subset inside the plotting script
        input_tau="../outputs/gi_scores/screen2023/construct_scores/all_gis_Tau.OI.Avg.tsv"
    output:
        output_figure_1e="../outputs/manuscript_figures/figure_1e.svg"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_1e.log"
    script:
        "../scripts/manuscript_figures/plot_figure_1e.R"

rule plot_figure_2a:
    input:
        input_nu_r1="../outputs/gi_scores/screen2023/differential_scores/gene_differential_scores_OI.R1.tsv",
        input_nu_r2="../outputs/gi_scores/screen2023/differential_scores/gene_differential_scores_OI.R2.tsv",
        input_nu_gene_level="../outputs/gi_scores/screen2023/differential_scores/differential_hits_OI.Avg.tsv"
    output:
        output_figure_2a="../outputs/manuscript_figures/figure_2a.png",
        output_figure_2a_labels="../outputs/manuscript_figures/figure_2a_labels.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_2a.log"
    script:
        "../scripts/manuscript_figures/plot_figure_2a.R"

rule plot_figure_2b:
    input:
        input_nu_r1="../outputs/gi_scores/screen2023/differential_scores/gene_differential_scores_OI.R1.tsv",
        input_nu_r2="../outputs/gi_scores/screen2023/differential_scores/gene_differential_scores_OI.R2.tsv",
        input_nu_gene_level="../outputs/gi_scores/screen2023/differential_scores/differential_hits_OI.Avg.tsv",
        input_idmap="data/annotations/screen2023_id_to_name_mapping.tsv"
    output:
        output_figure_2b="../outputs/manuscript_figures/figure_2b.png",
        output_figure_2b_labels="../outputs/manuscript_figures/figure_2b_labels.png"
    log:
        "../outputs/logs/manuscript_figures/plot_figure_2b.log"
    script:
        "../scripts/manuscript_figures/plot_figure_2b.R"


# rule cluster_genetic_interaction_scores:
#     input:
#         input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
#     output:
#         output_clustered_scores="../outputs/gi_scores/{screen}/clustered_scores/clustered_gene_combination_scores_{score}.tsv",
#         output_dendrogram="../outputs/gi_scores/{screen}/clustered_scores/dendrogram_{score}.rds"
#     log:
#         "../outputs/logs/{screen}/{screen}_{score}_cluster_genetic_interaction_scores.log"
#     params:
#         score=lambda wildcards: wildcards.score,
#         screen=lambda wildcards: wildcards.screen
#     script:
#         "../scripts/cluster_gi_scores.R"

# rule create_circos_plots:
#     input:
#         input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
#         input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
#     output:
#         output_circos_plot="../outputs/gi_scores/{screen}/plots/circos_plot_{score}.png"
#     log:
#         "../outputs/logs/{screen}/{screen}_{score}_create_circos_plots.log"
#     params:
#         score=lambda wildcards: wildcards.score,
#         screen=lambda wildcards: wildcards.screen
#     script:
#         "../scripts/create_circos_plots.R"