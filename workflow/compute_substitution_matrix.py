# Retrieve information from Snakemake
# We have to use check if values are strings here because inputs and outputs 
# are always iterables when they are defined with the syntax
# rules.<rule>.<input/output>, even if there is only one input / output.
input_file_path = snakemake.input if isinstance(snakemake.input, str) else snakemake.input[0]
output_file_path = snakemake.output if isinstance(snakemake.output, str) else snakemake.output[0]

# Initialize a dictionary where results will be stored
# Structure will be a dictionary of dictionaries:
# substitution_counts --> {ref_allele: {alt_allele: count}}
substitution_counts = {}

# Iterate over variant lines in VCF file and parse them into the
# results dictionary
with open(input_file_path) as vcf_file:
    for line in vcf_file:
        if line.startswith("#"):  # Ignore comments and header
            continue
        fields = line.rstrip('\n').split('\t')  # Split variant lines by tab
        ref = fields[3]  # Ref allele is in the 4th field (0-based)
        alt = fields[4]  # Alt allele is in the 3rd field (0-based)
        # We only count variants if ref and alt are different
        # and we only process SNPs (single nucleotide variants)
        if ref != alt and len(ref) == 1 and len(alt) == 1:
            if ref not in substitution_counts:
                # If it's the first time we encounter this reference nucleotide,
                # we create a new empty dictionary as value for this key
                substitution_counts[ref] = {}
            if alt not in substitution_counts[ref]:
                # If it's the first time we encounter this alt nucleotide for the
                # given ref nucleotide, we initialize the counter to 0
                substitution_counts[ref][alt] = 0
            # Increment counter for this combination of ref and alt alleles
            substitution_counts[ref][alt] += 1

# Iterate over results dictionary to save the results in comma-separated table
# We sort ref nucleotides and alt nucleotides so the order of rows and columns is the same
with open(output_file_path, 'w') as substitution_matrix:
    for ref_allele, alt_alleles in sorted(substitution_counts.items()):
        # We use the 'sep'.join() syntax to create a string of 'sep'-separated counts
        substitution_matrix.write(','.join([str(count) for nucleotide, count in sorted(alt_alleles.items())]) + '\n')
