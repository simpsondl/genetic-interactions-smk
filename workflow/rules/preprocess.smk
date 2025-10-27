rule apply_count_filters:
    input:
        input_counts="data/counts/{screen}_raw_counts.tsv"
    output:
        output_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    params:
        individual_sgRNA_median_threshold=config["INDIVIDUAL_SGRNA_MEDIAN_THRESHOLD"],
        combination_sgRNA_count_threshold=config["COMBINATION_SGRNA_COUNT_THRESHOLD"],
        counts_cols=lambda wildcards: config[f"{wildcards.screen.upper()}_COUNTS_COLUMNS"]
    script:
        "../scripts/apply_filters.R"

rule calculate_phenotypes:
    input:
        input_counts="data/counts/{screen}_raw_counts.tsv",
        input_filter_flags="../outputs/misc_results/{screen}_filter_flags.tsv"
    output:
        output_phenotypes="../outputs/phenotypes/{screen}_phenotypes.tsv",
        output_orientation_indep_phenotypes="../outputs/phenotypes/{screen}_orientation_independent_phenotypes.tsv",
        output_single_sgRNA_phenotypes="../outputs/phenotypes/{screen}_single_sgRNA_phenotypes.tsv"
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
    params:
        no_correlation_threshold=config["NO_CORRELATION_THRESHOLD"]
    script:
        "../scripts/apply_correlation_filter.R"