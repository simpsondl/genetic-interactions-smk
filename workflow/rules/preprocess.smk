rule apply_count_filters:
    input:
        input_counts=_choose_counts
    output:
        output_filter_flags="../outputs/misc_results/{screen,[^_]+}_filter_flags.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_apply_count_filters.log"
    params:
        individual_sgRNA_median_threshold=config["INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD"],
        combination_sgRNA_count_threshold=config["COMBINATION_SGRNA_COUNT_THRESHOLD"],
        counts_cols=lambda wildcards: config[f"{wildcards.screen.upper()}_COUNTS_COLUMNS"]
    script:
        "../scripts/apply_filters.R"

rule calculate_phenotypes:
    input:
        input_counts=_choose_counts,
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        output_phenotypes="../outputs/phenotypes/{screen,[^_]+}_phenotypes.tsv",
        output_orientation_indep_phenotypes="../outputs/phenotypes/{screen,[^_]+}_orientation_independent_phenotypes.tsv",
        output_single_sgRNA_phenotypes="../outputs/phenotypes/{screen,[^_]+}_single_sgRNA_phenotypes.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_calculate_phenotypes.log"
    params:
        counts_cols=lambda wildcards: config[f"{wildcards.screen.upper()}_COUNTS_COLUMNS"],
        pseudocount=config["PSEUDOCOUNT"],
        normalize=config["NORMALIZE"],
        doublings=lambda wildcards: config[f"{wildcards.screen.upper()}_DOUBLINGS"]
    script:
        "../scripts/calculate_phenotypes.R"

rule apply_correlation_filter:
    input:
        input_phenotypes="../outputs/phenotypes/{screen}_phenotypes.tsv",
        input_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_orientation_independent_phenotypes.tsv",
        input_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_single_sgRNA_phenotypes.tsv",
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        output_filtered_phenotypes="../outputs/phenotypes/{screen}_filtered_phenotypes.tsv",
        output_filtered_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_filtered_single_sgRNA_phenotypes.tsv",
        output_full_filter_flags="../outputs/misc_results/{screen}_full_filter_flags.tsv",
        output_correlation_results="../outputs/misc_results/{screen}_correlation_results.tsv"
    log:
        "../outputs/logs/{screen}/{screen}_apply_correlation_filter.log"
    params:
        no_correlation_threshold=config["NO_CORRELATION_THRESHOLD"]
    script:
        "../scripts/apply_correlation_filter.R"