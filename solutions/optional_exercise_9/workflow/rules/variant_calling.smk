# Implement rules to call variants from a set of BAM files and compute a matrix of
# substitution rates from the resulting VCF.


rule call_variants:
    '''
    Call variants from alignment files for all samples defined in the list SAMPLES.
    Sorted BAMs for all samples are collected with expand().
    The value of the substitution rate parameter -P in bcftools call is defined in
    params.
    '''
    input:
        alignment_files = expand(rules.samtools_sort.output, sample=config['data']['samples']),
        assembly = 'data/genome.fa'
    output:
        'results/variants.vcf'
    params:
        substitution_rate = config['bcftools']['mutation_rate']
    log:
        'logs/call_variants.txt'
    conda:
        # Here, the path to the environment file is relative to the CURRENT FILE
        '../envs/variant_calling_env.yaml'
    shell:
        'bcftools mpileup -f {input.assembly} {input.alignment_files} 2> {log}| '
        'bcftools call -P {params.substitution_rate} -mv - > {output} 2>> {log}'


rule compute_substitution_matrix:
    '''
    Compute a matrix of substitution frequencies from the VCF file generated by call_variants.
    The matrix is computed using a Python script that we call with Snakemake.
    '''
    input:
        rules.call_variants.output
    output:
        'results/substitution_matrix.tsv'
    conda:
        # Here, the path to the environment file is relative to the CURRENT FILE
        '../envs/variant_calling_env.yaml'  # Allows to specify the exact version of Python to use
    script:
        # Here, the path to the script file is relative to the CURRENT FILE
        '../scripts/compute_substitution_matrix.py'