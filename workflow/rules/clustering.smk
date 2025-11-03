rule diagnostic_plot:
    input:
        input_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_combination_scores_{{score}}.tsv",
        input_idmap=lambda wildcards: (
            f"data/annotations/{wildcards.screen}_id_to_name_mapping.tsv" 
            if config.get("COUNTS_SOURCE") == "manuscript" 
            else f"{OUTPUTS_DIR}/annotations/{wildcards.screen}_gene_combination_id_map.tsv"
        )
    output:
        output_diagnostic_plot=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/clusters/diagnostic_plot_{{score}}.svg"
    log:
        f"{LOGS_DIR}/{{screen}}/clustering/{{screen}}_{{score}}_diagnostic_plot.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/diagnostic_plot.R"

###########################################
# Wrapper rules to perform all clustering #
###########################################

rule plot_all_diagnostic_plots:
    input:
        lambda wildcards: _expand_diagnostic_plots(wildcards)