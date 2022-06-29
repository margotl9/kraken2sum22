resources = "/Users/margotlavitt/resources/"
results = "/Users/margotlavitt/Sum22/kraken2/results/"


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


# Align reads to virus catalog using bowtie2
rule build_bowtie2_db:
    input:
        # results
        # + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
        resources + "mgv_db/kraken_formatted_mgv.fasta"
    output:
        resources + "bowtie2_db/virusdb.1.bt2",
    params:
        db=resources + "bowtie2_db/virusdb",
    conda:
        "workflow/envs/bowtie2.yml"
    # threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads 1
        """
 
 
# Align reads to virus catalog using bowtie2
rule bowtie2:
    input:
        R1=resources + "input/kraken_test_reads/enriched_ERR2155200_R1.fastq.gz",
        R2=resources + "input/kraken_test_reads/enriched_ERR2155200_R2.fastq.gz",
        # R1=results
        # + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_1.fastq",
        # R2=results
        # + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_2.fastq",
        # R1S=results
        # + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_1.fastq",
        # R2S=results
        # + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_2.fastq",
        db=resources + "bowtie2_db/virusdb.1.bt2",
    output:
        results + "bowtie2_align/bam_files/bowtie2.bam",
    params:
        db=resources + "bowtie2_db/virusdb",
        sam=results + "bowtie2_align/bowtie2.sam",
    conda:
        "workflow/envs/bowtie2.yml"
    # threads: config["virus_abundance"]["metapop_threads"]
    # -U {input.R1S},{input.R2S} \
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads 1 \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {params.sam}\
 
        # convert sam to bam
        samtools view -S -b {params.sam} > {output}
        rm {params.sam}
        """


# storing kraken rules
#
# # build kraken database including custom virus database
# rule kraken_build:
#     input:
#         resources + "mgv_db/kraken_formatted_mgv.fasta",
#     output:
#         resources + "mgv_kraken2db/hash.k2d",
#     params:
#         db=resources + "mgv_kraken2db/",
#     conda:
#         "workflow/envs/kraken2.yml"
#     shell:
#         """
#         kraken2-build --download-taxonomy --db {params.db}
#         kraken2-build --add-to-library {input} --db {params.db}
#         kraken2-build --build --db {params.db}
#         """
#
# #align reads to kraken database
# rule kraken2:
#     input:
#         db=resources + "mgv_kraken2db/hash.k2d",
#         R1=resources + "input/kraken_test_reads/enriched_ERR2155200_R1.fastq.gz",  #forward reads
#         R2=resources + "input/kraken_test_reads/enriched_ERR2155200_R2.fastq.gz",  #reverse seq reads
#     output:
#         classification="results/kraken2.kraken",
#         report="results/kraken2.kreport",
#     params:
#         db=resources + "mgv_kraken2db/",
#     conda:
#         "workflow/envs/kraken2.yml"
#     shell:
#         """
#         kraken2 --paired {input.R1} {input.R2} \
#         --db {params.db} \
#         --report {output.report} > {output.classification}
#         """


rule bracken_build:
    input:
        resources + "mgv_kraken2db/hash.k2d",
    output:
        resources + "mgv_kraken2db/database150mers.kmer_distrib",
    params:
        db=resources + "mgv_kraken2db",
    threads: 1
    conda:
        "workflow/envs/bracken.yml"
    shell:
        """
        bracken-build -d {params.db} -t {threads} -k 35 -l 150
        """


rule bracken:
    input:
        db=resources + "mgv_kraken2db/database150mers.kmer_distrib",
        report="results/kraken2.kreport",
    output:
        "results/bracken_abundances.bracken",
    params:
        db=resources + "mgv_kraken2db",
    conda:
        "workflow/envs/bracken.yml"
    shell:
        """
        bracken -d {params.db} \
        -i {input.report} \
        -o {output} \
        -r 150 -l F -t 10
        """


rule all:
    input:
        "results/bracken_abundances.bracken",
