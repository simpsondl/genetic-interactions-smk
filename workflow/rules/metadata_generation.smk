rule generate_id_maps:
    input:
        input_counts="data/counts/{screen}_counts_only.tsv"
    output:
        output_construct_map="../outputs/annotations/{screen}_construct_id_map.tsv",
        output_guidecombination_map="../outputs/annotations/{screen}_guide_combination_id_map.tsv",
        output_genecombination_map="../outputs/annotations/{screen}_gene_combination_id_map.tsv",
        output_pseudogene_map="../outputs/annotations/{screen}_pseudogene_id_map.tsv",
        output_counts_w_metadata="../outputs/counts/{screen}_counts_with_metadata.tsv"
    log:
        "../outputs/logs/{screen}/metadata_generation/{screen}_generate_id_maps.log"
    conda:
        "../envs/smk-env.yaml"
    script:
        "../scripts/generate_id_maps.R"

rule validate_counts:
    input:
        input_counts=lambda wildcards: _choose_counts(wildcards)
    output:
        output_validation_marker=temp("../outputs/misc_results/{screen}_counts_validated.txt")
    log:
        "../outputs/logs/{screen}/metadata_generation/{screen}_validate_counts.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        expected_count_columns=config.get("EXPECTED_COUNT_COLUMNS", None)
    script:
        "../scripts/validate_counts_with_metadata.R"

