#!/usr/bin/env python3

from collections import OrderedDict

import argparse

parser = argparse.ArgumentParser(description='Concatenate several RSEM matrices')
parser.add_argument("-i", "--input-files", help="RSEM count input files", nargs='+', required=True)
parser.add_argument("-s", "--statistic", help="Column to extract (expected_count, TPM or FPKM)",
                    default="expected_count", choices=["expected_count", "TPM", "FPKM"])
parser.add_argument("-o", "--output", help="Output file name", default="rsem_counts.tsv")

args = parser.parse_args()
output_file = args.output  # Output file (concat matrices)
iofiles = args.input_files  # List of input rsem count files (matrices)
rsem_statistic = args.statistic  # Column on the matrices to extract for each file

nb_files = len(iofiles)
result_lines = OrderedDict()
with open(output_file, "w") as rsem_summary_file:
    # Header :
    rsem_summary_file.write("gene_id\ttranscript_id(s)")
    for i in range(0, nb_files):  # For each input file
        rs_file = iofiles.pop(0)
        print("Load file %d/%d: %s" % (i + 1, nb_files, rs_file))
        rsem_summary_file.write("\t" + ("count_%d" % i))  # Add file to header
        # Get count_content
        with open(rs_file, "r") as rsem_file:
            lines = rsem_file.readlines()
            header = lines[0].strip("\n").split("\t")
            tpm_col = header.index(rsem_statistic)
            for j in range(1, len(lines)):
                line = lines[j].strip("\n").split("\t")
                if (line[0], line[1]) not in result_lines:
                    result_lines[(line[0], line[1])] = []
                if len(result_lines[(line[0], line[1])]) < i:
                    for k in range(0, i):
                        result_lines[(line[0], line[1])].append("")
                result_lines[(line[0], line[1])].append(line[tpm_col])

    rsem_summary_file.write("\n")
    for h_line, v_line in result_lines.items():
        rsem_summary_file.write("".join([h_line[0], "\t", h_line[1], "\t", "\t".join(v_line), "\n"]))
