rule call_variants:
    '''
    Generate a pileup file from all the alignment files and directly
    call variants from this pileup.
    '''
    input:
        reference = config['reference'],
        alignments = expand('results/{sample}.sorted.bam', sample=config["samples"]),
        indexes = expand('results/{sample}.sorted.bam.bai', sample=config["samples"])
    output:
        'results/variants.vcf'
    params:
        subrate = config['subrate_bcftools']
    log:
        'logs/call_variants.log'
    benchmark:
        'benchmarks/call_variants.txt'
    conda:
        '../envs/variants.yaml'
    shell:
        'bcftools mpileup -f {input.reference} {input.alignments} | '
        'bcftools call -P {params.subrate} -mv - > {output} 2>{log}'


rule compute_substitution_table:
    '''
    Compute a substitution table from the VCF output of call_variants.
    '''
    input:
        rules.call_variants.output
    output:
        'results/substitution_table.tsv'
    script:
        '../scripts/create_substitution_table.py'
