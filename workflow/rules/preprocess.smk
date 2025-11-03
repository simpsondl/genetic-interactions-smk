rule apply_count_filters:
    input:
        input_counts=f"{OUTPUTS_DIR}/counts/{{screen}}_counts_with_metadata.tsv",
        validated_counts=f"{OUTPUTS_DIR}/misc_results/{{screen}}_counts_validated.txt"
    output:
        output_filter_flags=f"{OUTPUTS_DIR}/misc_results/{{screen,[^_]+}}_filter_flags.tsv"
    log:
        f"{LOGS_DIR}/{{screen}}/preprocess/{{screen}}_apply_count_filters.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        individual_sgRNA_median_threshold=lambda wildcards: config.get(f"{wildcards.screen.upper()}_INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD",
                                                                        config.get("INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD")),
        combination_sgRNA_count_threshold=lambda wildcards: config.get(f"{wildcards.screen.upper()}_COMBINATION_SGRNA_COUNT_THRESHOLD",
                                                                        config.get("COMBINATION_SGRNA_COUNT_THRESHOLD")),
        # counts_prefixes and start_of_screen may be provided per-screen or globally
        counts_prefixes=lambda wildcards: config.get(f"{wildcards.screen.upper()}_COUNT_PREFIXES", 
                                                     config.get("COUNT_PREFIXES", None)),
        start_of_screen=lambda wildcards: config.get(f"{wildcards.screen.upper()}_INITIAL_CONDITION_ID", 
                                                     config.get("INITIAL_CONDITION_ID", None))
    script:
        "../scripts/apply_filters.R"

rule calculate_phenotypes:
    input:
        input_counts=f"{OUTPUTS_DIR}/counts/{{screen}}_counts_with_metadata.tsv",
        input_filter_flags=f"{OUTPUTS_DIR}/misc_results/{{screen}}_filter_flags.tsv"
    output:
        # intermediate raw phenotype file; marked temp so it is removed after downstream rules
        output_phenotypes_raw=temp(f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_phenotypes.raw.tsv")
    log:
        f"{LOGS_DIR}/{{screen}}/preprocess/{{screen}}_calculate_phenotypes.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        pseudocount=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PSEUDOCOUNT",
                                                  config.get("PSEUDOCOUNT")),
        # doublings, replicates, phenotype_prefix_map, counts_prefixes may be provided per-screen or globally
        doublings=lambda wildcards: config[f"{wildcards.screen.upper()}_DOUBLINGS"],
        # replicates, phenotype_prefix_map, counts_prefixes may be provided per-screen or globally
        replicates=lambda wildcards: config.get(f"{wildcards.screen.upper()}_REPLICATES", 
                                                config.get("REPLICATES", None)),
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", 
                                                          config.get("PHENOTYPE_PREFIX_MAP", None)),
        counts_prefixes=lambda wildcards: config.get(f"{wildcards.screen.upper()}_COUNT_PREFIXES", 
                                                     config.get("COUNT_PREFIXES", None)),
        
    script:
        "../scripts/calculate_phenotypes.R"

rule normalize_phenotypes:
    input:
        input_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen}}_phenotypes.raw.tsv"
    output:
        # intermediate normalized phenotype file; still temporary
        output_normalized_phenotypes=temp(f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_phenotypes.normalized.tsv")
    log:
        f"{LOGS_DIR}/{{screen}}/preprocess/{{screen}}_normalize_phenotypes.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        normalize=lambda wildcards: config.get(f"{wildcards.screen.upper()}_NORMALIZE",
                                                config.get("NORMALIZE")),
        doublings=lambda wildcards: config[f"{wildcards.screen.upper()}_DOUBLINGS"],
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", config.get("PHENOTYPE_PREFIX_MAP", None))
    script:
        "../scripts/normalize_phenotypes.R"

rule create_averaged_phenotypes:
    input:
        input_phenotypes_normalized=f"{OUTPUTS_DIR}/phenotypes/{{screen}}_phenotypes.normalized.tsv",
        input_filter_flags=f"{OUTPUTS_DIR}/misc_results/{{screen}}_filter_flags.tsv"
    output:
        output_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_processed_phenotypes.tsv",
        output_orientation_indep_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_orientation_independent_phenotypes.tsv",
        output_single_sgRNA_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_single_sgRNA_phenotypes.tsv"
    log:
        f"{LOGS_DIR}/{{screen}}/preprocess/{{screen}}_create_averaged_phenotypes.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        average_replicates=lambda wildcards: config.get(f"{wildcards.screen.upper()}_AVERAGE_REPLICATES",
                                                        config.get("AVERAGE_REPLICATES")),
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", config.get("PHENOTYPE_PREFIX_MAP", None))
    script:
        "../scripts/create_averaged_phenotypes.R"

# Get phenotypes list for correlation filtering (only those that need filtering)
# Use NO_CORRELATION_FILTER config to determine which phenotypes get filtered output files
FILTER_PHENOS = config.get("NO_CORRELATION_FILTER", [])

rule apply_correlation_filter:
    input:
        input_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen}}_processed_phenotypes.tsv",
        input_orientation_indep_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen}}_orientation_independent_phenotypes.tsv",
        input_single_sgRNA_phenotypes=f"{OUTPUTS_DIR}/phenotypes/{{screen}}_single_sgRNA_phenotypes.tsv",
        input_filter_flags=f"{OUTPUTS_DIR}/misc_results/{{screen}}_filter_flags.tsv"
    output:
        # Dynamically create per-phenotype outputs only for phenotypes in NO_CORRELATION_FILTER
        **{f"output_filtered_{p.lower()}_phenotypes": f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_filtered_{p.lower()}_phenotypes.tsv" for p in FILTER_PHENOS},
        **{f"output_filtered_{p.lower()}_single_sgRNA_phenotypes": f"{OUTPUTS_DIR}/phenotypes/{{screen,[^_]+}}_filtered_{p.lower()}_single_sgRNA_phenotypes.tsv" for p in FILTER_PHENOS},
        output_full_filter_flags=f"{OUTPUTS_DIR}/misc_results/{{screen,[^_]+}}_full_filter_flags.tsv",
        output_correlation_results=f"{OUTPUTS_DIR}/misc_results/{{screen,[^_]+}}_correlation_results.tsv",
        output_correlation_summary=f"{OUTPUTS_DIR}/misc_results/{{screen,[^_]+}}_correlation_summary.tsv"
    log:
        f"{LOGS_DIR}/{{screen,[^_]+}}/preprocess/apply_correlation_filter.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        # which phenotypes to process (list) - can be screen-specific or global
        phenotypes=lambda wildcards: config.get(f"{wildcards.screen.upper()}_GI_PHENOTYPES",
                                                config.get("GI_PHENOTYPES", ["Gamma","Tau"])),
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", config.get("PHENOTYPE_PREFIX_MAP", None)),
        replicates=lambda wildcards: config.get(f"{wildcards.screen.upper()}_REPLICATES", config.get("REPLICATES", None)),
        # mapping of phenotype -> boolean indicating whether to run no-correlation filtering for that phenotype
        no_correlation_filter=lambda wildcards: config.get(f"{wildcards.screen.upper()}_NO_CORRELATION_FILTER",
                                                            config.get("NO_CORRELATION_FILTER", {})),
        no_correlation_threshold=lambda wildcards: config.get(f"{wildcards.screen.upper()}_NO_CORRELATION_THRESHOLD",
                                                               config.get("NO_CORRELATION_THRESHOLD")),
        correlation_filter_mode=lambda wildcards: config.get(f"{wildcards.screen.upper()}_CORRELATION_FILTER_MODE",
                                                              config.get("CORRELATION_FILTER_MODE", "avg_only"))
    script:
        "../scripts/apply_correlation_filter.R"