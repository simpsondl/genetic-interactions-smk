rule generate_id_maps:
    input:
        input_counts="data/counts/{screen}_counts_only.tsv"
    output:
        output_construct_map=f"{OUTPUTS_DIR}/annotations/{{screen}}_construct_id_map.tsv",
        output_guidecombination_map=f"{OUTPUTS_DIR}/annotations/{{screen}}_guide_combination_id_map.tsv",
        output_genecombination_map=f"{OUTPUTS_DIR}/annotations/{{screen}}_gene_combination_id_map.tsv",
        output_pseudogene_map=f"{OUTPUTS_DIR}/annotations/{{screen}}_pseudogene_id_map.tsv",
        output_counts_w_metadata=f"{OUTPUTS_DIR}/counts/{{screen}}_counts_with_metadata.tsv"
    log:
        f"{LOGS_DIR}/{{screen}}/metadata_generation/{{screen}}_generate_id_maps.log"
    conda:
        "../envs/smk-env.yaml"
    script:
        "../scripts/generate_id_maps.R"

rule validate_counts:
    input:
        input_counts=lambda wildcards: _choose_counts(wildcards)
    output:
        output_validation_marker=temp(f"{OUTPUTS_DIR}/misc_results/{{screen}}_counts_validated.txt")
    log:
        f"{LOGS_DIR}/{{screen}}/metadata_generation/{{screen}}_validate_counts.log"
    conda:
        "../envs/smk-env.yaml"
    params:
        expected_meta_columns=config.get("EXPECTED_META_COLUMNS", None)
    script:
        "../scripts/validate_counts_with_metadata.R"

