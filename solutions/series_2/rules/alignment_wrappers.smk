def get_sample_path(wildcards):
    '''
    Simple one-line function to return the path to a sample file from
    the value of the wildcard 'sample'.
    '''
    return config['samples'][wildcards.sample]


rule align:
    '''
    Align reads to an assembly using the BWA mem Snakemake wrapper
    '''
    input:
        reads = get_sample_path
    output:
        temp('results/{sample}.bam')
    threads:
        config['thread_bwa']
    log:
        'logs/{sample}_align.log'
    benchmark:
        'benchmarks/{sample}_align.txt'
    params:
        index = config['reference']
    wrapper:
        '0.78.0/bio/bwa/mem'


rule sort_bam:
    '''
    Sort the alignment files generate by the 'align' rule using the
    Samtools sort Snakemake wrapper.
    '''
    input:
        rules.align.output
    output:
        'results/{sample}.sorted.bam'
    log:
        'logs/{sample}_sort_bam.log'
    benchmark:
        'benchmarks/{sample}_sort_bam.txt'
    wrapper:
        '0.78.0/bio/samtools/sort'


rule index_bam:
    '''
    Index the sorted alignment files from the 'sort_bam' rule using the
    Samtools index Snakemake wrapper.
    '''
    input:
        rules.sort_bam.output
    output:
        'results/{sample}.sorted.bam.bai'
    log:
        'logs/{sample}_index_bam.log'
    benchmark:
        'benchmarks/{sample}_index_bam.txt'
    wrapper:
        '0.78.0/bio/samtools/index'
