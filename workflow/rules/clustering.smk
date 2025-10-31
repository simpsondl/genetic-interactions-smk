rule diagnostic_plot:
    input:
        input_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_diagnostic_plot="../outputs/gi_scores/{screen}/clusters/diagnostic_plot_{score}.svg"
    log:
        "../outputs/logs/{screen}/{screen}_{score}_diagnostic_plot.log"
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