#!/usr/bin/env python2

from Bio import SeqIO   # Not used in current state, but possibly will be required in future iterations

from collections import defaultdict


# Argument: SeqIO-parsed GenBank record
# Return value: list of bases, each populated with a list of gene(s) it is a part of
def build_bases_list_from_gb_record(record):

    # Ensure that the 0th element of Features is 'source', as per GenBank convention        
    if record.features[0].type == 'source':
        # Initialize an empty list with the same amount of elements as the total
        # number of base pairs, which is the length specified by the source
        bases = [[] for i in range(record.features[0].location.end)]
    else:
        raise AttributeError("Error: 0th element of Bio.SeqFeature.SeqFeature is not 'source' as expected")

    # Build bases[], which will store which gene each base belongs to, by iterating over all features
    # within the GenBank file (source, gene, CDS, etc) (~10,000 of these for contig_1.gb)
    for feature in record.features:

        # Iterate over CDS and RNA features because they contain gene information and more. 
        # Note: 'regulatory' features do not have a locus tag.
        if feature.type in ['CDS','tRNA','tmRNA','rRNA','ncRNA','regulatory']:

            # Iterate over all base indices for this gene (about 5,000,000 of these for contig_1.gb)
            for base_index in range(
                    feature.location.start+2,
                    feature.location.end-1):
                
                # Add 'features' dictionary to the corresponding base pair index 
                # Note that more than one gene can be stored in a given index of 
                # bases[], which allows us to account for overlapping genes
                # Also, indices have been adjusted to shrink "gene" to only where Himar1 can disrupt.
                bases[base_index].append(feature)
    
    return bases



# Argument: list of bases, each populated with a list of gene(s) it is a part of, such as the return of build_bases_list()
# Return value: list of overlapping genes lists
def build_overlapping_genes_list(bases):

    overlapping_genes_list = []

    # Iterate over all bases
    for base_gene_list in bases:
        
        # If the current base has more than one gene, it has overlapping genes
        if len(base_gene_list) > 1:

            # Create a list of overlapping genes, which will be added to a list of lists below if unique
            overlapping_genes = sorted(base_gene_list)

            # Create a new list of lists of overlapping genes, which is a list of lists of features
            # ("gene dictionaries"), and if the list above is not in our new list, add it (preserve uniqueness)
            if overlapping_genes not in overlapping_genes_list:
                overlapping_genes_list.append(overlapping_genes)
    
    return overlapping_genes_list



# Argument: list of overlapping genes lists, such as the return of build_overlapping_genes_list()
# Return value: dictionary keyed by the amount of overlaps, which maps to the amount of occurrences
def build_overlap_stats_dict(overlapping_genes_list):
    
    num_overlaps_stats = defaultdict(int)
   
    for overlapping_genes in overlapping_genes_list:
        # Count the overlapping genes by initializing the dictionary
        # when appropriate and calculate the number of overlaps
        num_overlaps = overlapping_genes[0].location.end - overlapping_genes[1].location.start
        num_overlaps_stats[num_overlaps] += 1

    return num_overlaps_stats
