"""
This part of the workflow handles transforming the data into standardized
formats and expects input file

    sequences_ndjson = "data/sequences.ndjson"

This will produce output files as

    metadata = "data/metadata.tsv"
    sequences = "data/sequences.fasta"

Parameters are expected to be defined in `config.transform`.
"""


rule fetch_general_geolocation_rules:
    output:
        general_geolocation_rules="data/general-geolocation-rules.tsv",
    params:
        geolocation_rules_url=config["transform"]["geolocation_rules_url"],
    shell:
        """
        # (1) Pick curl or wget based on availability    
        if which curl > /dev/null; then
            download_cmd="curl -fsSL --output"
        elif which wget > /dev/null; then
            download_cmd="wget -O"
        else
            echo "ERROR: Neither curl nor wget found. Please install one of them."
            exit 1
        fi

        # (2) Fetch general geolocation rules
        $download_cmd {output.general_geolocation_rules} {params.geolocation_rules_url}
        """


rule concat_geolocation_rules:
    input:
        general_geolocation_rules="data/general-geolocation-rules.tsv",
        local_geolocation_rules=config["transform"]["local_geolocation_rules"],
    output:
        all_geolocation_rules="data/all-geolocation-rules.tsv",
    shell:
        """
        cat {input.general_geolocation_rules} {input.local_geolocation_rules} >> {output.all_geolocation_rules}
        """


rule transform:
    input:
        sequences_ndjson="data/sequences_{serotype}.ndjson",
        all_geolocation_rules="data/all-geolocation-rules.tsv",
    output:
        metadata="data/raw_metadata_{serotype}.tsv",
        sequences="data/sequences_{serotype}.fasta",
    log:
        "logs/transform_{serotype}.txt",
    params:
        field_map=config["transform"]["field_map"],
        strain_regex=config["transform"]["strain_regex"],
        strain_backup_fields=config["transform"]["strain_backup_fields"],
        date_fields=config["transform"]["date_fields"],
        expected_date_formats=config["transform"]["expected_date_formats"],
        articles=config["transform"]["titlecase"]["articles"],
        abbreviations=config["transform"]["titlecase"]["abbreviations"],
        titlecase_fields=config["transform"]["titlecase"]["fields"],
        authors_field=config["transform"]["authors_field"],
        authors_default_value=config["transform"]["authors_default_value"],
        abbr_authors_field=config["transform"]["abbr_authors_field"],
        annotations=config["transform"]["annotations"],
        annotations_id=config["transform"]["annotations_id"],
        metadata_columns=config["transform"]["metadata_columns"],
        id_field=config["transform"]["id_field"],
        sequence_field=config["transform"]["sequence_field"],
    shell:
        """
        # (1) Pick curl or wget based on availability    
        if which curl > /dev/null; then
            download_cmd="curl -fsSL --output"
        elif which wget > /dev/null; then
            download_cmd="wget -O"
        else
            echo "ERROR: Neither curl nor wget found. Please install one of them."
            exit 1
        fi

        # (2) Download the required scripts if not already present
        [[ -d bin ]] || mkdir bin
        [[ -f bin/transform-field-names ]]      || $download_cmd bin/transform-field-names "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-field-names"
        [[ -f bin/transform-string-fields ]]    || $download_cmd bin/transform-string-fields "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-string-fields"
        [[ -f bin/transform-strain-names ]]     || $download_cmd bin/transform-strain-names "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-strain-names"
        [[ -f bin/transform-date-fields ]]      || $download_cmd bin/transform-date-fields "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-date-fields"
        [[ -f bin/transform-genbank-location ]] || $download_cmd bin/transform-genbank-location "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-genbank-location"
        [[ -f bin/transform-authors ]]          || $download_cmd bin/transform-authors "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/transform-authors"
        [[ -f bin/apply-geolocation-rules ]]    || $download_cmd bin/apply-geolocation-rules "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/apply-geolocation-rules"
        [[ -f bin/merge-user-metadata ]]        || $download_cmd bin/merge-user-metadata "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/merge-user-metadata"
        [[ -f bin/ndjson-to-tsv-and-fasta ]]    || $download_cmd bin/ndjson-to-tsv-and-fasta "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/ndjson-to-tsv-and-fasta"
        chmod +x bin/*

        # (3) Transform the data
        (cat {input.sequences_ndjson} \
            | ./bin/transform-field-names \
                --field-map {params.field_map} \
            | ./bin/transform-string-fields --normalize \
            | ./bin/transform-strain-names \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | ./bin/transform-date-fields \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | ./bin/transform-genbank-location \
            | ./bin/transform-string-fields \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | ./bin/transform-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | ./bin/apply-geolocation-rules \
                --geolocation-rules {input.all_geolocation_rules} \
            | ./bin/merge-user-metadata \
                --annotations {params.annotations} \
                --id-field {params.annotations_id} \
            | ./bin/ndjson-to-tsv-and-fasta \
                --metadata-columns {params.metadata_columns} \
                --metadata {output.metadata} \
                --fasta {output.sequences} \
                --id-field {params.id_field} \
                --sequence-field {params.sequence_field} ) 2>> {log}
        """

rule post_process_metadata:
    input:
        metadata="data/raw_metadata_{serotype}.tsv",
    output:
        metadata="data/metadata_{serotype}.tsv",
    shell:
       """
       ./bin/post_process_metadata.py --metadata {input.metadata} --outfile {output.metadata}
       """

rule compress:
    input:
        sequences="data/sequences_{serotype}.fasta",
        metadata="data/metadata_{serotype}.tsv",
    output:
        sequences="data/sequences_{serotype}.fasta.zst",
        metadata="data/metadata_{serotype}.tsv.zst",
    shell:
        """
        zstd -T0 -o {output.sequences} {input.sequences}
        zstd -T0 -o {output.metadata} {input.metadata}
        """