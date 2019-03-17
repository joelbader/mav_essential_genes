#!/usr/bin/env python2

import os, sys
from collections import OrderedDict

def enumerate_stat_names(path):
    
    if os.path.isdir(path):
        # Use first .tn_stats file it comes across
        for f in sorted(os.listdir(path)):
            if f.endswith(".tn_stats"):
                tn_stats_filename = os.path.join(path, f)
                break
    elif os.path.isfile(path):
        tn_stats_filename = path
    else:
        raise FileNotFoundError
    
    stat_names = []
    
    f = open(tn_stats_filename, "r")
    for line in f:
        
        if not line.startswith("# "):
            continue
        if line.startswith("#  "):
            continue

        stripped_line = line[2:].strip().split(':')[0].split(' ')[0]
        stat_names.append(stripped_line)

    return stat_names



def get_stats(path, num_contigs, stat_names):
    
    if os.path.isfile(path):
        single_file_mode_filename = os.path.basename(path)
        path = os.path.dirname(path)
    elif os.path.isdir(path):
        single_file_mode_filename = None
    else:
        raise FileNotFoundError

    stats = OrderedDict()
    
    for stat_name in stat_names:

        stats[stat_name] = OrderedDict()
        
        # Iterate over all samples' .tn_stats files
        for filename in sorted(os.listdir(path)):
            
            if single_file_mode_filename is not None and filename != single_file_mode_filename:
                continue
                
            if not filename.endswith(".tn_stats"):
                continue
                
            path_to_file = os.path.join(path, filename)
            num_lines_read_since_stat = 0
            sample_name = os.path.splitext(os.path.basename(filename))[0]

            f = open(path_to_file, "r")
            for line in f:
                
                if not line.startswith("# "):
                    continue
                
                stripped_line = line[2:].strip()
                if stat_name in stripped_line:
                    
                    # If the line doesn't end with a colon, the statistic is reported on the same line, and is not unique between the replicons
                    if not stripped_line.endswith(':'):
                        if ": " in stripped_line:
                            entry = stripped_line.split(": ")[-1]
                        else:
                            entry = stripped_line.split(' ')[1:]
                        stats[stat_name][sample_name] = []
                        stats[stat_name][sample_name].append(entry)
                        break
                    else:
                        stats[stat_name][sample_name] = []
                        num_lines_read_since_stat += 1
                elif num_lines_read_since_stat > 0:
                    num_lines_read_since_stat += 1
                    stats[stat_name][sample_name].append(line.split()[-1])
                if num_lines_read_since_stat is (num_contigs + 1):
                    break
            f.close()

    return stats


def main():
    if len(sys.argv) > 1:
        if sys.argv[1][0] == '/':
            path = sys.argv[1]
        else:
            path = os.path.join(os.getcwd(), sys.argv[1])
        
        valid_stat_names = enumerate_stat_names(path)
        
        if len(sys.argv) > 3:
            num_contigs = int(sys.argv[2])
            stat_names = sys.argv[3:]
            for stat_name in stat_names:
                if stat_name not in valid_stat_names:
                    print("Error: invalid stat name %s specified" % stat_name)
                    print("Available statistics:")
                    for name in valid_stat_names:
                        print("  %s" % name)
                    return
            stats = get_stats(path, num_contigs, stat_names)
            for stat_name in stats:
                if len(stats) > 1:
                    print("===== %s =====" % stat_name)
                for sample in stats[stat_name]:
                    print("%s: %s" % (sample, '\t '.join(stats[stat_name][sample])))
                if len(stats) > 1:
                    print("")
        else:
            print("Available statistics:")
            for name in valid_stat_names:
                print("  %s" % name)
    else:
        print("Usage: %s <dir_of_tn_stats_files> <num_replicons> [stat_name] [stat_name]" % sys.argv[0])
        print("  Example: %s 3 ./ TAs_hit template_count" % sys.argv[0])
        print("To enumerate the possible statistic names, run it with just the path argument:")
        print("  Example: %s ./tpp_output" % sys.argv[0])


if __name__ == "__main__":
    main()
