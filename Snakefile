# conda: "envs/nextstrain.yaml"
configfile: "config/config_zika.yaml"

rule all:
    input:
        auspice_json = "auspice/zika.json",

rule files:
    params:
        dropped_strains = "config/dropped_strains.txt",
        include_strains = "config/include_strains.txt",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",
        description = "config/description.md"

files = rules.files.params

rule download:
    message: "Downloading sequences and metadata from data.nextstrain.org"
    output:
        sequences = "data/sequences.fasta.zst",
        metadata = "data/metadata.tsv.zst"
    params:
        sequences_url = "https://data.nextstrain.org/files/zika/sequences.fasta.zst",
        metadata_url = "https://data.nextstrain.org/files/zika/metadata.tsv.zst"
    shell:
        """
        curl -fsSL --compressed {params.sequences_url:q} --output {output.sequences}
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        """

rule decompress:
    message: "Decompressing sequences and metadata"
    input:
        sequences = "data/sequences.fasta.zst",
        metadata = "data/metadata.tsv.zst"
    output:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv"
    shell:
        """
        zstd -d -c {input.sequences} > {output.sequences}
        zstd -d -c {input.metadata} > {output.metadata}
        """

rule append_local_sequences:
    message:
        """
        Combine and deduplicate aligned FASTAs from multiple origins.
        """
    input:
        decompressed_remote_input = rules.decompress.output.sequences,
        additional_local_inputs = lambda w: [origin["sequences"] for origin in config.get("inputs", {})]
    output:
        sequences = "results/combined_sequences.fasta"
    log:
        "logs/combine_input_sequences.txt"
    shell:
        """
        cat {input.additional_local_inputs} {input.decompressed_remote_input} > {output.sequences} 
        """


rule append_local_metadata:
    message:
        """
        Appending additional local input metadata.
        """
    input:
        decompressed_remote_input = rules.decompress.output.metadata
    output:
        metadata = "results/combined_metadata.tsv"
    log:
        "logs/combine_input_metadata.txt"
    run:
        import pandas as pd
        origins = ["nextstrain"] + [origin["name"] for origin in config.get("inputs", {})]
        metadatas = [input.decompressed_remote_input] + [origin["metadata"] for origin in config.get("inputs", {})]
        merged_metadata = pd.DataFrame()
        for origin, metadata in zip(origins, metadatas):
            metadata_df = pd.read_csv(metadata, sep='\t').assign(origin=origin)
            merged_metadata = pd.concat([merged_metadata, metadata_df])
        merged_metadata.to_csv(output.metadata, na_rep='?', sep='\t', index=False)


rule wrangle_metadata:
    input:
        metadata="data/metadata.tsv" if "inputs" not in config else "results/combined_metadata.tsv",
    output:
        metadata="results/wrangled_metadata.tsv",
    shell:
        """
        csvtk -t rename -f strain -n strain_original {input.metadata} \
          | csvtk -t mutate -f accession -n strain > {output.metadata}
        """


rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length} (50% of Zika virus genome)
        """
    input:
        sequences = rules.decompress.output.sequences if "inputs" not in config else "results/combined_sequences.fasta",
        metadata = rules.wrangle_metadata.output.metadata,
        exclude = files.dropped_strains,
        include = files.include_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country year month",
        sequences_per_group = 40,
        min_date = 2012,
        min_length = 5385
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
    output:
        alignment = "results/aligned.fasta"
    params:
        ref = config["reference"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-name {params.ref} \
            --output {output.alignment} \
            --fill-gaps \
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.wrangle_metadata.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
    output:
        node_data = "results/aa_muts.json",
        alignments = expand("results/translations/{gene}.fasta", gene=["ENV"])
    params:
        reference_annotation = config["reference_annotation"]
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --genes ENV \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {params.reference_annotation} \
            --output-node-data {output.node_data} \
            --alignment-output results/translations/%GENE.fasta
        """

rule additive_phenotype_prediction:
    input:
        alignment = "results/translations/ENV.fasta"
    output:
        node_data = "results/{experiment}_variant_phenotype_prediction.json",
        pred_data = "results/{experiment}_variant_phenotype_prediction.csv"
    log:
        "logs/{experiment}_additive_phenotype_prediction.txt"
    params:
        dms_wt_seq_id = lambda w: config["diff_sel_models"][f"{w.experiment}"]["dms_wt_seq_id"],
        mut_effects_df = lambda w: config["diff_sel_models"][f"{w.experiment}"]["mut_effects_df"],
        mut_effect_col = lambda w: config["diff_sel_models"][f"{w.experiment}"]["mut_effect_col"],
        mutation_col = lambda w: config["diff_sel_models"][f"{w.experiment}"]["mutation_col"],
    conda:
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        python my_profiles/dmsa-pred/dmsa_pred.py phenotype-prediction \
            --model-type additive \
            --alignment {input.alignment} \
            --dms-wt-seq-id {params.dms_wt_seq_id} \
            --mut-effects-df {params.mut_effects_df} \
            --mut-effect-col {params.mut_effect_col} \
            --mutation-col {params.mutation_col} \
            --experiment-label {wildcards.experiment} \
            --output-json {output.node_data} \
            --output-df {output.pred_data} 2>&1 | tee {log}
        """


def _get_dms_prediction_data(wildcards):

    inputs = []
    
    if "diff_sel_models" in config:
        inputs += list(expand(
            rules.additive_phenotype_prediction.output.node_data,
            experiment=list(config["diff_sel_models"])
        ))

    return inputs

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.wrangle_metadata.output.metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "region country",
        sampling_bias_correction = 3
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """



rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.wrangle_metadata.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        phenotype_predictions = _get_dms_prediction_data, 
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config,
        description = files.description
    output:
        auspice_json="auspice/zika.json",
        root_sequence="auspice/zika_root-sequence.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.phenotype_predictions} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "data ",
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
