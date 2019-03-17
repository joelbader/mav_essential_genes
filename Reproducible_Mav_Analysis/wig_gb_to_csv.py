#!/usr/bin/env python2

from Bio import SeqIO
from bio_parser_tools import build_bases_list_from_gb_record

import os, time, sys, argparse
from collections import OrderedDict

# Argument: entry of a CSV file to be escaped/sanitized
# Return value: escaped/sanitized version of the argument string
def csv_escape(csv_entry):
    # Remove all non-printable characters, note this does not remove non-printable unicode characters
    csv_entry = ''.join([c for c in csv_entry if ord(c) > 31 or ord(c) == 9])
    
    # Replace double-quotes with double double-quotes
    csv_entry.replace('"', '""')

    # If a comma is present in the csv entry, enclose the csv entry in double quotes
    if ',' in csv_entry:
        csv_entry = "\"%s\"" % (csv_entry)

    return csv_entry



# Helper function for create_csv_from_wig_and_genes(), to prevent code reuse, as it is run twice:
#     once for the unique identifier, and over all fields for the other one
# Arguments: bases list of list of genes, the base index, which fields to enumerate, the current  
#            fields columns for the CSV, and whether to run in unique mode (stop on first field found)
# Return value: list of fields' values for the given base_index
def get_gene_fields_list(bases, base_index, fields, current_fields_columns, record_id, unique_mode=False):
   
    # If fields list is empty, just return the current fields list
    if not fields:
        return current_fields_columns
    
    # Unique mode will cause early loop termination upon successful finding of a field
    if unique_mode:
        loop_break  = False
    
    # For every gene present at that base index, append the specified fields to the data
    for field in fields:
        if unique_mode:
            if loop_break:
                break

        gene_field_string = ""
        
        for gene in bases[base_index]:
            
            # Check if the gene has the specified field
            if field in gene.qualifiers:
                gene_field_string += gene.qualifiers[field][0] + "; "
                if unique_mode:
                    loop_break = True
            elif field == 'record_id_start_end_strand':
                if gene.location.strand is None:
                    strand_string = "N"
                elif gene.location.strand is 0:
                    strand_string = "U"
                elif gene.location.strand > 0:
                    strand_string = "+"
                elif gene.location.strand < 0:
                    strand_string = "-"
                record_id_start_end_strand_string = "%s_%d-%d_%s" % (record_id, gene.location.start + 1, gene.location.end, strand_string)
                gene_field_string += record_id_start_end_strand_string + "; "
                if unique_mode:
                    loop_break = True
            elif field == 'REF_start_end_strand':
                if gene.location.strand is None:
                    strand_string = "N"
                elif gene.location.strand is 0:
                    strand_string = "U"
                elif gene.location.strand > 0:
                    strand_string = "+"
                elif gene.location.strand < 0:
                    strand_string = "-"
                ref_start_end_strand_string = "REF_%d-%d_%s" % (gene.location.start + 1, gene.location.end, strand_string)
                gene_field_string += ref_start_end_strand_string + "; "
                if unique_mode:
                    loop_break = True
            elif field == 'strand':
                if gene.location.strand is None:
                    strand_string = "N/A"
                elif gene.location.strand is 0:
                    strand_string = "Unknown"
                elif gene.location.strand > 0:
                    strand_string = "+"
                elif gene.location.strand < 0:
                    strand_string = "-"
                gene_field_string += strand_string + "; "
            elif field == 'record_id':
                gene_field_string += record_id + "; "
        if unique_mode:
            if len(bases[base_index]) is 0:
                # Break out of the loop if there aren't any genes to iterate over
                loop_break = True

            if loop_break is True:
                # Append the column if it's the unique element has been found
                current_fields_columns.append(csv_escape(gene_field_string[:-2]))
        else:
            current_fields_columns.append(csv_escape(gene_field_string[:-2]))

    return current_fields_columns

# Argument: GenBank filename
# Return value: List containing the valid Record IDs from the GenBank file
def enumerate_record_ids_from_genbank(genbank_file):
    # Use BioPython's SeqIO to parse the GenBank file
    return sorted(SeqIO.to_dict(SeqIO.parse(genbank_file, "genbank")))

# Arguments: wig prefixes, genbank filename, record_id to specify, unique identifier, list of fields to add to file, and 
#            whether or not to include column names (should only be set to True once per CSV output)
# Return value: csv formatted buffer
def create_csv_from_wig_and_genes(wig_prefixes_list_filename, wig_prefixes, genbank_file, record_id, unique, fields, include_column_names=True, single_contig=False):
    
    if unique is None or not unique:
        unique = ['record_id_start_end_strand']
    
    # Remove duplicates from unique and fields
    args.unique = list(OrderedDict.fromkeys(args.unique))
    args.fields = list(OrderedDict.fromkeys(args.fields))

    # Enumerate files that satisfy the record id
    wig_files = []
    for prefix in wig_prefixes:
        
        # Identify path and true prefix
        if prefix[0] is '/':
            path = os.path.dirname(prefix)
        else:
            path = os.path.dirname(os.path.abspath(wig_prefixes_list_filename))
        
        # Identify suffix
        if single_contig:
            suffix = ".wig"
        else:
            suffix = "%s.wig" % record_id

        # Identify the file that satisfies both prefix and suffix and add it to list
        for filename in sorted(os.listdir(path)):
            if os.path.isfile(os.path.join(path, filename)) \
                    and os.path.basename(os.path.join(path, filename)).startswith(os.path.basename(prefix))  \
                    and filename.endswith(suffix):
                wig_files.append(os.path.join(path, filename))
                break
    if len(wig_files) is 0:
        raise FileNotFoundError

    with open(wig_files[0], "r") as wig:
        
        # Populate csv_line[] with list of each line, each line being a list from the file, delimited by the space character
        csv_line = []
        for line in wig:
            csv_line.append([record_id] + line.rstrip('\n').split(' '))
    
    
    # Use BioPython's SeqIO to parse the GenBank file
    genbank_dict = SeqIO.to_dict(SeqIO.parse(genbank_file, "genbank"))
    
        
    if (len(genbank_dict) is 1) and (record_id is None or record_id is ""):
        record_id = genbank_dict.keys()[0]
    
    if record_id not in genbank_dict:
        raise KeyError

    try:
        bases = build_bases_list_from_gb_record(genbank_dict[record_id])
    except AttributeError as err:
        raise

    # Iterate over all lines of the Wiggle file, with the exception of the first two (metadata)
    for wig_line_number in range(2, len(csv_line)):
        
        # Get the base index specified from the line from the Wiggle file
        base_index = int(csv_line[wig_line_number][1])
        
        # If base_index specified by Wiggle file is out of bounds, ignore that line and print warning.
        if base_index > len(bases) - 1:
            sys.stderr.write("[!] [WARNING]: Line %d of '%s' (line = '%s %s') specifies an out of bounds base index %d (max specified by the GenBank file is %d) for record %s. The line will be ignored and the CSV will continue to be created, but most likely this record is not meant for this Wiggle file..\n" 
                    % (wig_line_number, wig_files[0], csv_line[wig_line_number][1], csv_line[wig_line_number][2], base_index, len(bases) - 1, record_id))
        else:
            gene_unique_identifier = get_gene_fields_list(bases, base_index, unique, [], record_id, True)
            gene_fields_list = get_gene_fields_list(bases, base_index, fields, gene_unique_identifier, record_id, False)
            csv_line[wig_line_number].extend(gene_fields_list)

    # Build the CSV by joining the now-modified csv_line list
    csv_contents = ""
    
    # Make lists for storing the metadata of the Wiggle files
    metadata1 = []
    metadata2 = []
    metadata1.append(' '.join(csv_line[0]))
    metadata2.append(' '.join(csv_line[1]))

    if include_column_names:
        # Enumerate columns of this CSV file
        temp_csv_contents = ""  # Use of this as a temporary string storage space so we can get the metadata from the other Wiggle files as well
        if len(fields) > 0:
            temp_csv_contents += "contig,insertion_site,unique_identifier (%s),%s" % ( " or ".join(unique), ','.join(fields))
            for prefix in wig_prefixes:
                temp_csv_contents += ",read_count (%s)" % os.path.basename(prefix)
        else:
            temp_csv_contents += "contig,insertion_site,unique_identifier (%s)" % ( " or ".join(unique))
            for prefix in wig_prefixes:
                temp_csv_contents += ",read_count (%s)" % prefix

    for line_number in range(2, len(csv_line)):
        # Move second element (number of insertions/read_count) to last column
        csv_line[line_number].append(csv_line[line_number].pop(2))
    
    # Open the other Wiggle files and make new columns for their insertion/read counts
    for wig_file_index in range(1, len(wig_files)):
        with open(wig_files[wig_file_index], "r") as wig:
            # Parse the metadata lines
            metadata1.append(wig.next().rstrip('\n'))
            metadata2.append(wig.next().rstrip('\n'))
            
            # Add the insertion counts
            line_number = 2
            for line in wig:
                csv_line[line_number].append(line.rstrip('\n').split(' ')[1])
                line_number += 1
    
    # Comment out the lines appending to csv_contents to skip metadata sections
    #csv_contents += ','.join(metadata1) + '\n'
    if metadata2.count(metadata2[0]) == len(metadata2):
    #    csv_contents += metadata2[0] + '\n'
        pass
    else:
    #    csv_contents += ','.join(metadata2) + '\n'
        sys.stderr.write("[!] [WARNING]: Mismatches within the second line of metadata from the Wiggle files.\n") 

    if include_column_names:
        csv_contents += temp_csv_contents + '\n'

    for line_number in range(2, len(csv_line)):
        # Convert the list to CSV
        csv_contents += ','.join(csv_line[line_number]) + '\n'
    
    return csv_contents



def main(args):

    exit_flag = False

    if not os.path.isfile(args.gbk):
        print("[-] Error: '%s' does not exist." % (args.gbk))
        exit_flag = True
    
    if exit_flag:
        return
    
    valid_record_ids = enumerate_record_ids_from_genbank(args.gbk)
    
    # Ensure that there is only one record if no Record ID was specified
    if not args.record and len(valid_record_ids) > 1:
        print("[-] Error: No Record ID specified (-r flag), and multiple records exist in '%s'." % (args.gbk))
        print("....Valid Record IDs for '%s' are:" % (args.gbk))
        for record_id in valid_record_ids:
            print("....  %s" % (record_id))
        return

    print("[+] %s starting with parameters:" % (sys.argv[0]))
    print("....Wiggle prefix(es)         %s" % (', '.join(args.prefix)))
    print("....Record ID(s)              %s" % (', '.join(args.record)))
    print("....GenBank file:             %s" % (args.gbk))
    print("....Unique identifier fields: %s" % (', '.join(args.unique)))
    if args.fields:
        print("....Fields:                   %s" % (', '.join(args.fields)))
    print("....Output file:              %s" % (args.outfile))
    
    list_of_csv_files_contents = []
    try:
        start_time = time.time()
        columns_enumerated = False
        for record in args.record:
            print("[+] Combining Wiggles with prefix '%s' with Record ID '%s' from '%s'..." % (args.prefix, record, args.gbk))
            start_time_record = time.time()
            if columns_enumerated:
                csv_file_contents = create_csv_from_wig_and_genes(args.list, args.prefix, args.gbk, record, args.unique, args.fields, False, False) 
            else:
                csv_file_contents = create_csv_from_wig_and_genes(args.list, args.prefix, args.gbk, record, args.unique, args.fields, True, (len(args.record) == 1)) 
                columns_enumerated = True
            list_of_csv_files_contents.append(csv_file_contents) 
            end_time_record = time.time()
            print("[+] Finished processing Record ID '%s' successfully in %f seconds." 
                    % (record, end_time_record - start_time_record))
        end_time = time.time()
    except KeyError:
        print("[-] Specified Record ID '%s' is not valid for the Genbank file." % (args.record))
        print("....Valid Record IDs for '%s' are:" % (args.gbk))
        for record_id in valid_record_ids:
            print("....  %s" % (record_id))
        return
    
    print("[+] Finished processing all records successfully in %f seconds." % (end_time - start_time))

    print("[+] Writing CSV files...")
    with open(args.outfile, "w") as csv_file:
        for csv in list_of_csv_files_contents:
            csv_file.write(csv)
        print("[+] '%s' written." % (args.outfile))

    print("[+] Success. Exiting...")



def parse_args():
    parser = argparse.ArgumentParser(description="Create CSV file by combining Wiggle files with GenBank files")
    optional_args = parser._action_groups.pop()

    required_args = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional_args)

    required_args.add_argument("--prefix", "--pre", "--wp", "-p", type=str, 
            nargs='+', help="Specify Wiggle (.wig) file prefixes", required=False)
    required_args.add_argument("--gbk", "--gb", "--genbank", "-g", type=str, 
            help="Specify GenBank (.gb) file", required=True)
    required_args.add_argument("--fields", "--field", "--sections", "--section", "-f", "-s", type=str, 
            nargs='+', help="Specify sections to parse from GenBank file to put in .csv file. " \
                    "Sample fields include 'locus_tag', 'inference', 'note', 'codon_start', "   \
                    "'transl_table', 'product', 'protein_id', 'translation', and the custom "   \
                    "'strand', 'record_id_start_end_strand' and 'REF_start_end_strand'", required=False)
    
    optional_args.add_argument("--outfile", "--output", "--out", "-o", type=str, 
            help="Specify output file in csv format (default is wigs_combined_[gb-name].csv)", required=False)
    optional_args.add_argument("--unique", "--unique-identifier", "-u", type=str,
            nargs='+', help="Specify a field/section to use as the unique identifier. If more than one "    \
                    "is specified, they will be used in descending priority order. Default if none are "    \
                    "specified is 'record_id_start_end_strand'.", required=False)
    optional_args.add_argument("--record", "--record-id", "--rec", "-r", type=str, 
            nargs='+', help="Specify Record ID(s) to parse from the GenBank file. If this argument is "     \
                    "omitted, we assume all records from the GenBank file are to be used.", required=False)
    optional_args.add_argument("--list", "--list-of-prefixes", "--lwp", "-l", type=str, 
            help="Specify file to load line-separated list of Wiggle prefixes from (replaces --prefix)", required=False)

    args = parser.parse_args()

    # Specify default outfile name if none given to sigs_combined_[gb].csv
    if args.outfile is None:
        args.outfile = "wigs_combined_%s.csv" % os.path.basename(args.gbk) 

    if not args.record:
        ### We want to be able to not specify any records, so this code is deprecated
        ## If no record specified, check that only one record is in the GenBank file and use that
        #valid_record_ids = enumerate_record_ids_from_genbank(args.gbk)
        #if len(valid_record_ids) is 1:
        #    args.record = valid_record_ids
        args.record = enumerate_record_ids_from_genbank(args.gbk)

    # If the prefixes are specified from a file (--list), parse them as if they were specified from the --prefix argument
    if args.list:
        args.prefix = []
        path = os.path.dirname(args.list)  # used to apply relative path from where the list is located

        with open(args.list, "r") as prefixes:
            for line in prefixes:
                args.prefix.append(os.path.join(path, line.rstrip('\n')))
    
    # Ensure that we have 'record_id_start_end_strand' as a backup in any case
    default_unique = 'record_id_start_end_strand'
    if args.unique is None:
        args.unique = [default_unique]
    if args.fields is None:
        args.fields = []
    if args.unique and default_unique not in args.unique:
        args.unique.append(default_unique)
   

    return args



if __name__ == '__main__':
    args = parse_args()
    main(args)
