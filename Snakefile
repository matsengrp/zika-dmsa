if not config:
    configfile: "config/config_zika.yaml"

rule all:
    input:
        auspice_json = "auspice/zika.json",

rule files:
    params:
        input_fasta = "data/zika.fasta",
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/zika_reference.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",
        description = "config/description.md"

files = rules.files.params

rule download:
    message: "Downloading sequences and metadata from data.nextstrain.org"
    output:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
    params:
        sequences_url = "https://data.nextstrain.org/files/workflows/zika/test/sequences_all.fasta.zst",
        metadata_url = "https://data.nextstrain.org/files/workflows/zika/test/metadata_all.tsv.zst"
    shell:
        """
        curl -fsSL {params.sequences_url:q} --output {output.sequences}
        curl -fsSL {params.metadata_url:q} --output {output.metadata}
        """

rule wrangle_metadata:
    input:
        metadata=rules.download.output.metadata,
    output:
        metadata="results/wrangled_metadata.tsv",
    params:
        strain_id=lambda w: config.get("strain_id_field", "strain"),
        wrangle_metadata_url="https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/scripts/wrangle_metadata.py"
    shell:
        """
        if [[ ! -d bin ]]; then
          mkdir bin
        fi
        if [[ ! -f bin/wrangle_metadata.py ]]; then
          cd bin
          wget https://raw.githubusercontent.com/nextstrain/monkeypox/master/scripts/wrangle_metadata.py
          chmod 755 *
          cd ..
        fi

        python3 ./bin/wrangle_metadata.py --metadata {input.metadata} \
            --strain-id {params.strain_id} \
            --output {output.metadata}
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
        sequences = rules.download.output.sequences,
        metadata = rules.wrangle_metadata.output.metadata,
        exclude = files.dropped_strains
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
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
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
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

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
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config,
        description = files.description
    output:
        auspice_json = "results/raw_zika.json",
        root_sequence = "results/raw_zika_root-sequence.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

rule final_strain_name:
    input:
        auspice_json=rules.export.output.auspice_json,
        metadata=rules.wrangle_metadata.output.metadata,
        root_sequence=rules.export.output.root_sequence,
    output:
        auspice_json=rules.all.input.auspice_json,
        root_sequence="auspice/zika_root-sequence.json",
    params:
        display_strain_field=lambda w: config.get("display_strain_field", "strain"),
        set_final_strain_name_url="https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/scripts/set_final_strain_name.py",
    shell:
        """
        if [[ ! -d bin ]]; then
          mkdir bin
        fi
        cd bin
        [[ -f set_final_strain_name.py ]] || wget {params.set_final_strain_name_url}
        chmod 755 *
        cd ..

        python3 bin/set_final_strain_name.py \
            --metadata {input.metadata} \
            --input-auspice-json {input.auspice_json} \
            --display-strain-name {params.display_strain_field} \
            --output {output.auspice_json}

        cp {input.root_sequence} {output.root_sequence}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "data ",
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
