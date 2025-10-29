rule plot_models:
    input:
        input_gamma="../outputs/gi_scores/screen2023/individual_scores/Gamma.OI.Avg/PARP1_+_226595744.23-P1P2.txt",
        input_tau="../outputs/gi_scores/screen2023/individual_scores/Tau.OI.Avg/PARP1_+_226595744.23-P1P2.txt"
    output:
        output_gamma_plot="../outputs/gi_scores/screen2023/plots/sgrna_model_PARP1_+_226595744.23-P1P2_gamma.svg",
        output_tau_plot="../outputs/gi_scores/screen2023/plots/sgrna_model_PARP1_+_226595744.23-P1P2_tau.svg"
    log:
        "../outputs/logs/screen2023/screen2023_plot_models.log"
    script:
        "../scripts/plot_sgrna_model.R"


rule cluster_genetic_interaction_scores:
    input:
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    output:
        output_clustered_scores="../outputs/gi_scores/{screen}/clustered_scores/clustered_gene_combination_scores_{score}.tsv",
        output_dendrogram="../outputs/gi_scores/{screen}/clustered_scores/dendrogram_{score}.rds"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_cluster_genetic_interaction_scores.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/cluster_gi_scores.R"

rule create_circos_plots:
    input:
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_circos_plot="../outputs/gi_scores/{screen}/plots/circos_plot_{score}.png"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_create_circos_plots.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/create_circos_plots.R"