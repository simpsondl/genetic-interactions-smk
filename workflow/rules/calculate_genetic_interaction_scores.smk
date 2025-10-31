rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes=_gi_orientation_indep_phenotype,
        input_single_sgRNA_phenotypes=_gi_single_sgRNA_phenotype
    output:
        output_dir=directory("../outputs/gi_scores/{screen}/individual_scores/{score}"),
        output_all_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        output_model_estimates="../outputs/gi_scores/{screen}/models/model_estimates_{score}.tsv",
        output_model_stats="../outputs/gi_scores/{screen}/models/model_stats_{score}.tsv",
        output_workspace=temp("../outputs/gi_scores/{screen}/construct_scores/gi_workspace_{score}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score}_compute_genetic_interaction_scores.log"
    conda:
        "../envs/smk-env.yaml"
    wildcard_constraints:
        score="(Gamma\..*|Tau\..*)"
    params:
        screen=lambda wildcards: wildcards.screen,
        score=lambda wildcards: wildcards.score
    script:
        "../scripts/calculate_gi_scores.R"

rule calculate_gene_level_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_gi_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_{score}.rds",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        output_gene_level_workspace=temp("../outputs/gi_scores/{screen}/gene_combination_scores/gene_level_workspace_{score}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score}_calculate_gene_level_scores.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_gene_level_scores.R"

rule calculate_discriminant_scores:
    input:
        input_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_{score}.tsv",
        input_gene_level_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_{score}.tsv",
        input_gene_level_workspace="../outputs/gi_scores/{screen}/gene_combination_scores/gene_level_workspace_{score}.rds"
    output:
        output_discriminant_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv",
        output_discriminant_workspace=temp("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_workspace_{score}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score}_calculate_discriminant_scores.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_discriminant_scores.R"

rule calculate_differential_scores:
    input:
        input_gamma_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Gamma.{rep}.tsv",
        input_tau_gi_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Tau.{rep}.tsv",
        input_gamma_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_Gamma.{rep}.rds",
        input_tau_workspace="../outputs/gi_scores/{screen}/construct_scores/gi_workspace_Tau.{rep}.rds",
        input_idmap="data/annotations/{screen}_id_to_name_mapping.tsv"
    output:
        output_differential_scores="../outputs/gi_scores/{screen}/construct_scores/all_gis_Nu.{rep}.tsv",
        output_gene_differential_scores="../outputs/gi_scores/{screen}/gene_combination_scores/gene_combination_scores_Nu.{rep}.tsv",
        output_discriminant_differential_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_Nu.{rep}.tsv",
        output_diff_workspace=temp("../outputs/gi_scores/{screen}/discriminant_scores/discriminant_workspace_Nu.{rep}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_Nu.{rep}_calculate_differential_scores.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_differential_scores.R"

rule call_hits:
    input:
        input_scores="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_scores_{score}.tsv",
        input_discriminant_workspace="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_workspace_{score}.rds"
    output:
        output_hits="../outputs/gi_scores/{screen}/discriminant_scores/discriminant_hits_{score}.tsv",
        output_hits_workspace=temp("../outputs/gi_scores/{screen}/discriminant_scores/hits_workspace_{score}.rds")
    log:
        "../outputs/logs/{screen}/{screen}_{score}_call_hits.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        threshold=lambda wildcards: config["DIFFERENTIAL_HIT_THRESHOLD"] if str(wildcards.score).startswith("Nu") else config["HIT_THRESHOLD"],
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/call_hits.R"

##############################################
# Wrapper rules to build all scores and hits #
##############################################

rule compute_all_genetic_interaction_scores:
    input:
        lambda wildcards: _expand_scores_for_screen(wildcards)

rule compute_all_gene_level_scores:
    input:
        lambda wildcards: _expand_gene_level_score_targets(wildcards)

rule compute_all_discriminant_scores:
    input:
        lambda wildcards: _expand_discriminant_score_targets(wildcards)

rule compute_all_differential_scores:
    input:
        lambda wildcards: _expand_differential_score_targets(wildcards)

rule identify_all_hits:
    input:
        lambda wildcards: _expand_hit_targets(wildcards)

# rule identify_differential_hits:
#     input:
#         lambda wildcards: _expand_differential_hit_targets(wildcards)
