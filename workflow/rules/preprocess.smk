rule apply_count_filters:
    input:
        input_counts=_choose_counts,
        validated_counts="../outputs/misc_results/{screen}_counts_validated.txt"
    output:
        output_filter_flags="../outputs/misc_results/{screen,[^_]+}_filter_flags.tsv"
    log:
        "../outputs/logs/{screen}/preprocess/{screen}_apply_count_filters.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        individual_sgRNA_median_threshold=config["INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD"],
        combination_sgRNA_count_threshold=config["COMBINATION_SGRNA_COUNT_THRESHOLD"],
        # counts_prefixes and start_of_screen may be provided per-screen or globally
        counts_prefixes=lambda wildcards: config.get(f"{wildcards.screen.upper()}_COUNT_PREFIXES", 
                                                     config.get("COUNT_PREFIXES", None)),
        start_of_screen=lambda wildcards: config.get(f"{wildcards.screen.upper()}_INITIAL_CONDITION_ID", 
                                                     config.get("INITIAL_CONDITION_ID", None))
    script:
        "../scripts/apply_filters.R"

rule calculate_phenotypes:
    input:
        input_counts=_choose_counts,
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        # intermediate raw phenotype file; marked temp so it is removed after downstream rules
        output_phenotypes_raw=temp("../outputs/phenotypes/{screen,[^_]+}_phenotypes.raw.tsv")
    log:
        "../outputs/logs/{screen}/preprocess/{screen}_calculate_phenotypes.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        pseudocount=config["PSEUDOCOUNT"],
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
        input_phenotypes="../outputs/phenotypes/{screen,[^_]+}_phenotypes.raw.tsv"
    output:
        # intermediate normalized phenotype file; still temporary
        output_normalized_phenotypes=temp("../outputs/phenotypes/{screen,[^_]+}_phenotypes.normalized.tsv")
    log:
        "../outputs/logs/{screen}/preprocess/{screen}_normalize_phenotypes.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        normalize=config["NORMALIZE"],
        doublings=lambda wildcards: config[f"{wildcards.screen.upper()}_DOUBLINGS"],
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", config.get("PHENOTYPE_PREFIX_MAP", None))
    script:
        "../scripts/normalize_phenotypes.R"

rule replicate_averaging:
    input:
        input_phenotypes_normalized="../outputs/phenotypes/{screen,[^_]+}_phenotypes.normalized.tsv",
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        output_phenotypes="../outputs/phenotypes/{screen,[^_]+}_phenotypes.tsv",
        output_orientation_indep_phenotypes="../outputs/phenotypes/{screen,[^_]+}_orientation_independent_phenotypes.tsv",
        output_single_sgRNA_phenotypes="../outputs/phenotypes/{screen,[^_]+}_single_sgRNA_phenotypes.tsv"
    log:
        "../outputs/logs/{screen}/preprocess/{screen}_replicate_averaging.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        average_replicates=config["AVERAGE_REPLICATES"],
        phenotype_prefix_map=lambda wildcards: config.get(f"{wildcards.screen.upper()}_PHENOTYPE_PREFIX_MAP", config.get("PHENOTYPE_PREFIX_MAP", None))
    script:
        "../scripts/replicate_averaging.R"

rule apply_correlation_filter:
    input:
        input_phenotypes="../outputs/phenotypes/{screen}_phenotypes.tsv",
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_orientation_independent_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_single_sgRNA_phenotypes.tsv",
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        output_filtered_gamma_phenotypes="../outputs/phenotypes/{screen}_filtered_gamma_phenotypes.tsv",
        output_filtered_tau_phenotypes="../outputs/phenotypes/{screen}_filtered_tau_phenotypes.tsv",
        output_filtered_gamma_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_gamma_single_sgRNA_phenotypes.tsv",
        output_filtered_tau_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_tau_single_sgRNA_phenotypes.tsv",
        output_full_filter_flags="../outputs/misc_results/{screen}_full_filter_flags.tsv",
        output_correlation_results="../outputs/misc_results/{screen}_correlation_results.tsv"
    log:
        "../outputs/logs/{screen}/preprocess/{screen}_apply_correlation_filter.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        no_correlation_threshold=config["NO_CORRELATION_THRESHOLD"]
    script:
        "../scripts/apply_correlation_filter.R"