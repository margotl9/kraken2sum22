resources = "/Users/margotlavitt/resources/"


rule virome_db_download:
    output:
        genomes=resources + "mgv_db/mgv_votu_reps.fasta",
        metadata=resources + "mgv_db/mgv_metadata.tsv",
    shell:
        """
        # Download viral genome database (MGV)
        curl  https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.fna \
        -o {output.genomes}

        # Download viral genome metadata (MGV)
        curl https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/mgv_contig_info.tsv \
        -o  {output.metadata}
        """


rule customize_virus_headers:
    input:
        genomes=resources + "mgv_db/mgv_votu_reps.fasta",
        metadata=resources + "mgv_db/mgv_metadata.tsv",
    output:
        resources + "mgv_db/kraken_formatted_mgv.fasta",
    conda:
        "workflow/envs/jupyter.yml"
    notebook:
        "workflow/notebooks/customize_virus_headers.py.ipynb"


# build kraken database including custom virus database
rule kraken_build:
    input:
        resources + "mgv_db/kraken_formatted_mgv.fasta",
    output:
        resources + "mgv_kraken2db/hash.k2d",
    params:
        db=resources + "mgv_kraken2db/",
    conda:
        "workflow/envs/kraken2.yml"
    shell:
        """
        kraken2-build --download-taxonomy --db {params.db}
        kraken2-build --add-to-library {input} --db {params.db}
        kraken2-build --build --db {params.db}
        """


#align reads to kraken database
rule kraken2:
    input:
        db=resources + "mgv_kraken2db/hash.k2d",
        R1=resources + "input/reads_ncbi_genomes_R1.fastq.gz",  #forward reads
        R2=resources + "input/reads_ncbi_genomes_R2.fastq.gz",  #reverse seq reads
    output:
        classification="results/kraken2.kraken",
        report="results/kraken2.kreport",
    params:
        db=resources + "mgv_kraken2db/",
    conda:
        "workflow/envs/kraken2.yml"
    shell:
        """
        kraken2 --paired {input.R1} {input.R2} \
        --db {params.db} \
        --report {output.report} > {output.classification}
        """


# rule bracken_build:
#     input:
#         resources + "kraken2db/",
#         # resources + "kraken2db/database150mers.kmer_distrib"
#         #resources + "virusdb/combined.fasta" #create resources folder and define var
#     output:
#         resources + "bracken/virome_brackendb",
#     threads: 1
#     conda:
#         "envs/bracken.yml"
#     shell:
#         """
#         bracken-build --db {output} -t {threads} -k 35 -l 150
#         """


# rule bracken:
#     input:
#         db=resources + "kraken2db/",
#         report="results/kraken2.kreport",
#     output:
#         "results/bracken_abundances.bracken",
#     conda:
#         "workflow/envs/bracken.yml"
#     shell:
#         """
#         bracken -d {input.db} \
#         -i {input.report} \
#         -o {output}.bracken \
#         -r 150 -l P -t 10
#         """


rule all:
    input:
        "results/bracken_abundances.bracken",
