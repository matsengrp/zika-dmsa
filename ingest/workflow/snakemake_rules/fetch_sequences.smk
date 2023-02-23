"""
This part of the workflow handles fetching sequences from various sources.
Uses `config.sources` to determine which sequences to include in final output.

Currently only fetches sequences from GenBank, but other sources can be
defined in the config. If adding other sources, add a new rule upstream
of rule `fetch_all_sequences` to create the file `data/{source}.ndjson` or the
file must exist as a static file in the repo.

Produces final output as

    sequences_ndjson = "data/sequences.ndjson"

"""

def download_serotype(w):
    serotype = {
        'all': '64320',
    }
    return serotype[w.serotype]

rule fetch_from_genbank:
    output:
        genbank_ndjson="data/genbank_{serotype}.ndjson",
    params:
        ncbi_taxon_id=download_serotype,
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
        [[ -f bin/csv-to-ndjson ]]      || $download_cmd bin/csv-to-ndjson "https://raw.githubusercontent.com/nextstrain/monkeypox/644d07ebe3fa5ded64d27d0964064fb722797c5d/ingest/bin/csv-to-ndjson"
        [[ -f bin/genbank-url ]]        || $download_cmd bin/genbank-url "https://raw.githubusercontent.com/nextstrain/dengue/ca659008bfbe4b3f799e11ecd106a0b95977fe93/ingest/bin/genbank-url"
        [[ -f bin/fetch-from-genbank ]] || $download_cmd bin/fetch-from-genbank "https://raw.githubusercontent.com/nextstrain/dengue/ca659008bfbe4b3f799e11ecd106a0b95977fe93/ingest/bin/fetch-from-genbank"
        chmod +x bin/*
        
        # (3) Fetch the sequences
        ./bin/fetch-from-genbank {params.ncbi_taxon_id} > {output.genbank_ndjson}
        """


def _get_all_sources(wildcards):
    return [f"data/{source}_{wildcards.serotype}.ndjson" for source in config["sources"]]


rule fetch_all_sequences:
    input:
        all_sources=_get_all_sources,
    output:
        sequences_ndjson="data/sequences_{serotype}.ndjson",
    shell:
        """
        cat {input.all_sources} > {output.sequences_ndjson}
        """
