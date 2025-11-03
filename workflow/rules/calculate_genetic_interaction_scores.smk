rule compute_genetic_interaction_scores:
    input:
        input_orientation_indep_phenotypes=_gi_orientation_indep_phenotype,
        input_single_sgRNA_phenotypes=_gi_single_sgRNA_phenotype
    output:
        output_dir=directory(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/individual_scores/{{score}}"),
        output_all_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/all_gis_{{score}}.tsv",
        output_model_estimates=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/models/model_estimates_{{score}}.tsv",
        output_model_stats=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/models/model_stats_{{score}}.tsv",
        output_workspace=temp(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/gi_workspace_{{score}}.rds")
    log:
        f"{LOGS_DIR}/{{screen}}/calculate_gi_scores/{{screen}}_{{score}}_compute_genetic_interaction_scores.log"
    conda:
        "../envs/smk-env.yaml"
    wildcard_constraints:
        # Build regex dynamically from GI_PHENOTYPES config (e.g., "Gamma|Tau|Rho")
        score=r"(" + "|".join([rf"{p}\..*" for p in config.get("GI_PHENOTYPES", ["Gamma", "Tau"])]) + r")"
    params:
        score=lambda wildcards: wildcards.score
    script:
        "../scripts/calculate_gi_scores.R"

rule calculate_gene_level_scores:
    input:
        input_gi_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/all_gis_{{score}}.tsv",
        input_gi_workspace=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/gi_workspace_{{score}}.rds",
        input_idmap=f"{OUTPUTS_DIR}/annotations/{{screen}}_gene_combination_id_map.tsv"
    output:
        output_gene_level_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_combination_scores_{{score}}.tsv",
        output_gene_level_workspace=temp(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_level_workspace_{{score}}.rds")
    log:
        f"{LOGS_DIR}/{{screen}}/calculate_gi_scores/{{screen}}_{{score}}_calculate_gene_level_scores.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_gene_level_scores.R"

rule calculate_discriminant_scores:
    input:
        input_gi_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/all_gis_{{score}}.tsv",
        input_gene_level_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_combination_scores_{{score}}.tsv",
        input_gene_level_workspace=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_level_workspace_{{score}}.rds"
    output:
        output_discriminant_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_scores_{{score}}.tsv",
        output_discriminant_workspace=temp(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_workspace_{{score}}.rds")
    log:
        f"{LOGS_DIR}/{{screen}}/calculate_gi_scores/{{screen}}_{{score}}_calculate_discriminant_scores.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        score=lambda wildcards: wildcards.score,
        screen=lambda wildcards: wildcards.screen
    script:
        "../scripts/calculate_discriminant_scores.R"

# Define DIFFERENTIAL_COMPARISONS mapping and regex for wildcards
DIFFERENTIAL_COMPARISONS = config.get("DIFFERENTIAL_COMPARISONS", {"Nu": ["Tau", "Gamma"]})
_diff_keys = list(DIFFERENTIAL_COMPARISONS.keys())
if len(_diff_keys) == 0:
    _diff_keys = ["Nu"]
_comp_regex = "(" + "|".join(_diff_keys) + ")"

rule calculate_differential_scores:
    input:
        # reference/treated inputs are resolved from the DIFFERENTIAL_COMPARISONS mapping using wildcard 'comp'
        input_reference_gi_scores=lambda wildcards: f"{OUTPUTS_DIR}/gi_scores/{wildcards.screen}/construct_scores/all_gis_{DIFFERENTIAL_COMPARISONS[wildcards.comp][1]}.{wildcards.rep}.tsv",
        input_treated_gi_scores=lambda wildcards: f"{OUTPUTS_DIR}/gi_scores/{wildcards.screen}/construct_scores/all_gis_{DIFFERENTIAL_COMPARISONS[wildcards.comp][0]}.{wildcards.rep}.tsv",
        input_reference_workspace=lambda wildcards: f"{OUTPUTS_DIR}/gi_scores/{wildcards.screen}/construct_scores/gi_workspace_{DIFFERENTIAL_COMPARISONS[wildcards.comp][1]}.{wildcards.rep}.rds",
        input_treated_workspace=lambda wildcards: f"{OUTPUTS_DIR}/gi_scores/{wildcards.screen}/construct_scores/gi_workspace_{DIFFERENTIAL_COMPARISONS[wildcards.comp][0]}.{wildcards.rep}.rds",
        input_idmap=lambda wildcards: f"{OUTPUTS_DIR}/annotations/{wildcards.screen}_gene_combination_id_map.tsv"
    output:
        output_differential_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/construct_scores/all_gis_{{comp}}.{{rep}}.tsv",
        output_gene_differential_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/gene_combination_scores/gene_combination_scores_{{comp}}.{{rep}}.tsv",
        output_discriminant_differential_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_scores_{{comp}}.{{rep}}.tsv",
        output_diff_workspace=temp(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_workspace_{{comp}}.{{rep}}.rds")
    log:
        f"{LOGS_DIR}/{{screen}}/calculate_gi_scores/{{screen}}_{{comp}}.{{rep}}_calculate_differential_scores.log"
    conda:
        "../envs/smk-env.yaml"
    wildcard_constraints:
        comp=_comp_regex,
    params:
        rep=lambda wildcards: wildcards.rep,
        screen=lambda wildcards: wildcards.screen,
        reference_pheno=lambda wildcards: DIFFERENTIAL_COMPARISONS[wildcards.comp][1],
        treated_pheno=lambda wildcards: DIFFERENTIAL_COMPARISONS[wildcards.comp][0]
    script:
        "../scripts/calculate_differential_scores.R"

rule call_hits:
    input:
        input_scores=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_scores_{{score}}.tsv",
        input_discriminant_workspace=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_workspace_{{score}}.rds"
    output:
        output_hits=f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/discriminant_hits_{{score}}.tsv",
        output_hits_workspace=temp(f"{OUTPUTS_DIR}/gi_scores/{{screen}}/discriminant_scores/hits_workspace_{{score}}.rds")
    log:
        f"{LOGS_DIR}/{{screen}}/calculate_gi_scores/{{screen}}_{{score}}_call_hits.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        # Extract the score prefix (e.g., "Gamma", "Tau", "Nu" from "Gamma.OI.R1")
        # and look it up in HIT_THRESHOLDS for precise threshold mapping.
        # Support screen-specific HIT_THRESHOLDS and falls back to global version.
        threshold=lambda wildcards: config.get(f"{wildcards.screen.upper()}_HIT_THRESHOLDS", 
                                                config.get("HIT_THRESHOLDS", {})).get(
            str(wildcards.score).split(".")[0],
            config.get(f"{wildcards.screen.upper()}_HIT_THRESHOLD",
                       config.get("HIT_THRESHOLD", 0.995))
        ),
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
