resources = "/Users/margotlavitt/resources/"

# rule virome_db_download:
#     output:
#         resources + "virusdb/combined.fasta",
#     shell:
#         """
#         # Download viral genome database (MGV)
#         curl  https://portal.nersc.gov/MGV/MGV_v1.0_2021_07_08/mgv_votu_representatives.fna \
#         -o /Users/margotlavitt/resources/virusdb/combined.fasta
#         """


rule kraken2db_download:
    output:
        resources + "kraken2db/hash.k2d",
    params:
        k2db_targz=resources + "kraken2db/k2_viral_20220607.tar.gz",
        k2db_dir=resources + "kraken2db/",
    shell:
        """
        # Download viral minikraken database
        curl  https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220607.tar.gz \
        -o {params.k2db_targz}

        # Extract the tar.gz file
        tar -xzf /Users/margotlavitt/resources/kraken2db/k2_viral_20220607.tar.gz \
        -C {params.k2db_dir}
        """

# build kraken database including custom virus database
# rule kraken_build:
#     input:
#         resources + "virusdb/combined.fasta",
#     output:
#         resources + "kraken2/combined.k2d",
#     conda:
#         "envs/kraken2.yml"
#     shell:
#         """
#         kraken2-build --db combined
#         """

#align reads to kraken database
rule kraken2:
    input:
        db = resources + "kraken2db/hash.k2d",
        R1 = resources + "input/reads_ncbi_genomes_R1.fastq.gz", #forward reads
        R2 = resources + "input/reads_ncbi_genomes_R2.fastq.gz", #reverse seq reads
    output:
        classification = "results/kraken2.kraken",
        report = "results/kraken2.kreport",
    params:
        db = resources + "kraken2db/",
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
#         # resources + "kraken2db/database150mers.kmer_distrib"
#         #resources + "virusdb/combined.fasta" #create resources folder and define var
#     output:
#         resources + "bracken/virome_brackendb"
#     threads: 1
#     conda:
#         "envs/bracken.yml"
#     shell:
#         """
#         bracken-build --db {output} -t {threads} -k 35 -l 150
#         """


rule bracken:
    input:
        db = resources + "kraken2db/",
        report = "results/kraken2.kreport",
    output:
        "results/bracken_abundances.bracken",
    conda:
        "workflow/envs/bracken.yml"
    shell:
        """
        bracken -d {input.db} \
        -i {input.report} \
        -o {output}.bracken \
        -r 150 -l P -t 10
        """

rule all:
    input: 
        "results/bracken_abundances.bracken"
