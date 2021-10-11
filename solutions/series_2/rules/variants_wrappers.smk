rule mpileup:
    '''
    Generate a pileup file from all the alignment files using the BCF tools
    pileup Snakemake wrapper.
    '''
    input:
        ref = config['reference'],
        index = f'{config["reference"]}.fai',
        alignments = expand(rules.sort_bam.output, sample=config['samples']),
        indexes = expand(rules.index_bam.output, sample=config['samples'])
    output:
        pileup = 'results/pileup.bcf',
    log:
        'logs/mpileup.log',
    wrapper:
        '0.78.0/bio/bcftools/mpileup'


rule call_variants:
    '''
    Call variants from the pileup file generate by the 'mpileup' rule using
    the BCF call Snakemake wrapper.
    '''
    input:
        pileup = rules.mpileup.output
    output:
        calls = 'results/variants.vcf'
    params:
        caller = '-m',
        options = f'-P {config["subrate_bcftools"]} -v',
    log:
        'logs/call_variants.log'
    benchmark:
        'benchmarks/call_variants.txt'
    wrapper:
        '0.78.0/bio/bcftools/call'


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
