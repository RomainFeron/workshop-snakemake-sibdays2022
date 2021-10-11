def get_sample_path(wildcards):
    '''
    Simple one-line function to return the path to a sample file from
    the value of the wildcard 'sample'.
    '''
    return config['samples'][wildcards.sample]


rule align:
    '''
    Align reads to an assembly
    '''
    input:
        reference = config['reference'],
        samples = get_sample_path
    output:
        temp('results/{sample}.bam')
    threads:
        config['thread_bwa']
    log:
        'logs/{sample}_align.log'
    benchmark:
        'benchmarks/{sample}_align.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'bwa mem -t {threads} {input.reference} {input.samples} | '
        'samtools view -b > {output} 2>{log}'


rule sort_bam:
    '''
    Sort the alignment files generate by the 'align' rule.
    '''
    input:
        rules.align.output
    output:
        'results/{sample}.sorted.bam'
    log:
        'logs/{sample}_sort_bam.log'
    benchmark:
        'benchmarks/{sample}_sort_bam.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'samtools sort -O bam {input} > {output} 2>{log}'


rule index_bam:
    '''
    Index the sorted alignment files from the 'sort_bam' rule.
    '''
    input:
        rules.sort_bam.output
    output:
        'results/{sample}.sorted.bam.bai'
    log:
        'logs/{sample}_index_bam.log'
    benchmark:
        'benchmarks/{sample}_index_bam.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'samtools index {input} 2>{log}'
