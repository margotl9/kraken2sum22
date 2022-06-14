#build kraken database including custom virus database
rule kraken_build:
    input:
        resources + "virusdb/combined.fasta" #create resources folder and define var
    output:
        resources + "kraken2/combined.k2d"
    shell:
        """
        kraken2-build --db combined
        """
#align reads to kraken database
rule kraken2:
    input:  
        db = resources + "kraken2/combined.k2d",
        R1 = resources + "input/sequences_1.fastq", #seq reads
        R2 = resources + "input/sequences_2.fastq", #reverse seq reads
    output:
        "results/kraken2_classification.txt"
    shell:
        """
        kraken2 --paired {input.R1} {input.R2} --db {input.db}
        """

# rule braken_build:
#     input:
#         resources + "virusdb/combined.fasta" #create resources folder and define var
#     output:
#         resources + "braken/combined.k2d"
#     shell:
        
# rule braken:
#     input:
#         db = resources + "braken/combined.k2d",
#         kraken_class = "results/kraken2_classification.txt"
#     output:
#         "results/braken_abundances.txt"
#     shell:

rule all:
    input: 
        "results/braken_abundances.txt"