rule cluster_genetic_interaction_scores:
    input:
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv"
    output:
        output_clustered_scores="../outputs/gi_scores/{screen}/clustered_scores/clustered_gene_combination_scores_{score}.tsv",
        output_dendrogram="../outputs/gi_scores/{screen}/clustered_scores/dendrogram_{score}.rds"
    log:
        "../outputs/logs/{screen}_{score}_cluster_genetic_interaction_scores.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/cluster_gi_scores.R"

rule plot_sgrna_model:
    input:
        input_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score}.tsv"
    output:
        output_plot="../outputs/gi_scores/{screen}/plots/sgrna_model_{score}.png"
    log:
        "../outputs/logs/{screen}_{score}_plot_sgrna_model.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/plot_sgrna_model.R"

rule create_circos_plots:
    input:
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_circos_plot="../outputs/gi_scores/{screen}/plots/circos_plot_{score}.png"
    log:
        "../outputs/logs/{screen}_{score}_create_circos_plots.log"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/create_circos_plots.R"