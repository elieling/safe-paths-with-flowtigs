#!/usr/bin/python3

"""
Convert the output of the different validation tools into a LaTeX file.
"""

import sys, subprocess, hashlib, os, pathlib, json
import traceback
import re

hashdir = os.path.abspath(sys.argv[1])
genome_name_file_name = sys.argv[2]
graph_statistics_file_name = sys.argv[3]
bandage_png_name = sys.argv[4]
combined_eaxmax_plot_name = sys.argv[5]
output_file_name = sys.argv[6]
runtimes_file = sys.argv[7]
nodes_and_edges = sys.argv[8]
number_of_characters = sys.argv[9]

experiments = []

for i in range(10, len(sys.argv), 5):
    if i + 4 >= len(sys.argv):
        exit("Number of experiment parameters not divisible by 5")

    # shortname, omnitig algo file, QUAST dir, ContigValidator file, resources evaluation file
    experiments.append((sys.argv[i], sys.argv[i + 1], sys.argv[i + 2], sys.argv[i + 3], sys.argv[i + 4]))

def append_latex_table_second_column(table, appendix):
    try:
        if len(table) == 0:
            return appendix

        table_keys = set([line[:line.index("&")].strip() for line in table if '&' in line])
        appendix_keys = set([line[:line.index("&")].strip() for line in appendix if '&' in line])
        table_value_column_count = table[0].count("&")

        table_index = 0
        appendix_index = 0
        result = []

        def append_rows(row, appendix):
            row = row.strip()
            if row[-2:] == "\\\\":
                row = row[:-2] # Remove trailing new line backslashes
            appendix = appendix[appendix.index("&"):].strip()
            return row + appendix

        def new_row(appendix):
            return appendix[:appendix.index("&")] + (" & " * table_value_column_count) + appendix[appendix.index("&"):]

        def append_none(row):
            row = row.strip()
            if row[-2:] == "\\\\":
                row = row[:-2] # Remove trailing new line backslashes
            return row + " & N/A \\\\"

        while table_index < len(table) or appendix_index < len(appendix):
            if table_index < len(table):
                if '&' not in table[table_index]:
                    table_index += 1
                    continue
            if appendix_index < len(appendix):
                if '&' not in appendix[appendix_index]:
                    appendix_index += 1
                    continue

            if table_index < len(table) and appendix_index < len(appendix):
                table_line = table[table_index]
                table_key = table_line[:table_line.index("&")].strip()
                appendix_line = appendix[appendix_index]
                appendix_key = appendix_line[:appendix_line.index("&")].strip()

                if table_key == appendix_key:
                    result.append(append_rows(table_line, appendix_line))
                    table_index += 1
                    appendix_index += 1
                elif table_key in appendix_keys:
                    # Appendix contains something extra
                    result.append(new_row(appendix_line))
                    appendix_index += 1
                elif appendix_key in table_keys:
                    # Appendix misses something
                    result.append(append_none(table_line))
                    table_index += 1
                else:
                    row = appendix_key + (" & " * table_value_column_count) + "\\\\"
                    result.append(append_none(row))
                    appendix_index += 1
                    #sys.exit("Found completely mismatching keys: {} and {}".format(table_key, appendix_key))
            elif table_index < len(table):
                # Appendix misses something
                result.append(append_none(table[table_index]))
                table_index += 1
            elif appendix_index < len(appendix):
                # Appendix contains something extra
                result.append(new_row(appendix[appendix_index]))
                appendix_index += 1
            else:
                assert False

        return result
    except Exception:
        print("Cannot append to table\ntable:\n{}\nappendix:\n{}".format('\n'.join(table), '\n'.join(appendix)))

        raise

#########################
### Process name file ###
#########################

name_file = open(genome_name_file_name, 'r')
name_lines = name_file.readlines()
name_lines = [x.replace("_", "\\_").replace("{", "").replace("}", "") for x in name_lines]

###############################
### Process algorithm files ###
###############################

# def read_algorithm_file(prefix):
#   try:
#       algorithm_file = open(os.path.join(prefix, "compute_injectable_contigs.tex"))
#       algorithm_lines = algorithm_file.readlines()
#       return algorithm_lines
#   except:
#       print("Did not find algorithm file for '" + prefix + "'")
#       return []

headline = "Parameter"
algorithm_table = []
# for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
#   headline += " & " + label
#   table = read_algorithm_file(os.path.join(prefix, "injectable_contigs"))
#   algorithm_table = append_latex_table_second_column(algorithm_table, table)

algorithm_table = [headline + "\\\\ \\hline"] + algorithm_table


###########################
### Process QUAST files ###
###########################

def read_quast_file(filename):
    try:
        quast_file = open(filename, 'r')
        quast_lines = quast_file.readlines()
        quast_lines = [x.strip() for x in quast_lines]
        quast_lines = [x.replace("\\hline", "") for x in quast_lines]
        quast_lines = quast_lines[8:-4] # Remove LaTeX header and footer
        quast_lines = [line.replace("misassemblies caused by fragmented reference", "mis. caused by frag. ref.") for line in quast_lines]
        return quast_lines
    except:
        print("Did not find quast file for '" + filename + "'")
        return []

headline = "Parameter"
quast_table = []
for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
    headline += " & " + label
    table = read_quast_file(os.path.join(quast_directory, "report.tex"))
    quast_table = append_latex_table_second_column(quast_table, table)

quast_table = [headline + "\\\\ \\hline"] + quast_table

headline = "Parameter"
quast_misassemblies_table = []
for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
    headline += " & " + label
    table = read_quast_file(os.path.join(quast_directory, "contigs_reports", "misassemblies_report.tex"))
    quast_misassemblies_table = append_latex_table_second_column(quast_misassemblies_table, table)

quast_misassemblies_table = [headline + "\\\\ \\hline"] + quast_misassemblies_table


#####################################
### Process ContigValidator files ###
#####################################

# def read_contig_validator_file(prefix):
#   try:
#       contig_validator_file = open(prefix + ".contigvalidator", 'r')
#       contig_validator_lines = contig_validator_file.readlines()
#       contig_validator_lines = contig_validator_lines[1].split()[1:] # Remove header and file name column
#       contig_validator_lines[0] = "\\%exact & " + contig_validator_lines[0] + "\\%\\\\"
#       contig_validator_lines[1] = "\\%align & " + contig_validator_lines[1].replace("%", "\\%") + "\\\\"
#       contig_validator_lines[2] = "recall & " + contig_validator_lines[2].replace("%", "\\%") + "\\\\"
#       contig_validator_lines[3] = "precision & " + contig_validator_lines[3].replace("%", "\\%") + "\\\\"
#       return contig_validator_lines
#   except:
#       print("Did not find contig validator file for '" + prefix + "'")
#       return []

headline = "Parameter"
contig_validator_table = []
# for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
#   headline += " & " + label
#   table = read_contig_validator_file(prefix)
#   contig_validator_table = append_latex_table_second_column(contig_validator_table, table)

contig_validator_table = [headline + "\\\\ \\hline"] + contig_validator_table

##########################################
### Process CLI graph statistics files ###
##########################################

graph_statistics_table = []
if graph_statistics_file_name != "none": 
    try:
        graph_statistics_file = open(graph_statistics_file_name, 'r')
        graph_statistics_lines = graph_statistics_file.readlines()
        graph_statistics_table = ["Parameter & Value \\\\ \\hline"] + [x.strip() for x in graph_statistics_lines]
    except:
        print("Did not find graph statistics file '" + graph_statistics_file_name + "'")
        graph_statistics_table = []

#########################################
### Process resources evaluation file ###
#########################################

headline = "Parameter"
resources_table = []
# for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
#     headline += " & " + label
#     try:
#         with open(resources_evaluation_name, 'r') as resources_evaluation_file:
#             resources_evaluation = json.load(resources_evaluation_file)
#             time = resources_evaluation["total"]["time"]
#             mem = resources_evaluation["total"]["mem"] / (1024**2)
#             subtable = [f"time [s] & {time:.2f} \\\\", f"mem [GiB] & {mem:.2f} \\\\"]

#             for key, values in resources_evaluation.items():
#                 if key == "total":
#                     continue
#                 key = key.replace("_", "\\_")
#                 time = values["time"]
#                 mem = values["mem"] / (1024**2)
#                 subtable.append(f"{key} time & {time:.2f} \\\\")
#                 subtable.append(f"{key} mem & {mem:.2f} \\\\")

#             resources_table = append_latex_table_second_column(resources_table, subtable)
#     except:
#         print("Error processing resources evaluation file")
#         raise

resources_table = [headline + "\\\\ \\hline"] + resources_table

########################
### Build LaTeX file ###
########################

def table_header(caption, column_count):
    #header = """
    #\\begin{table}[ht]
    #\\begin{center}
    #\\caption{""" + caption + """}
    #\\begin{tabular}{|l*{1}{|r}|}
    #\\hline
    #"""
    header = """
    \\begin{table}[ht]
    \\begin{center}
    \\fontsize{6pt}{7pt}\\selectfont
    \\caption{\\fontsize{6pt}{7pt}\\selectfont\\ """ + caption + """}
    \\begin{tabular}{|l*{1}|"""
    for _ in range(column_count):
        header += "r"
    header += """|}
    \\hline
    """
    return header

table_footer = """\\hline
    \\end{tabular}
    \\end{center}
    \\end{table}
    """

def format_row_numbers(row):
    row = row.strip()
    row = row.rstrip('\\')
    result = []

    for column in row.split('&'):
        column = column.strip()
        if column.isdigit():
            column = int(column)
            column = f"{column:,}"
        else:
            try:
                column = float(column)
                column = f"{column:,.2f}"
            except ValueError:
                pass
        result.append(column)

    return " & ".join(result) + "\\\\"

def format_metatable_row_numbers(row):
    row = row.strip()
    row = row.rstrip('\\')
    result = []

    for column in row.split('&'):
        column = column.strip()
        if column.isdigit():
            column = int(column)
            column = f"{column:,}"
        else:
            try:
                column = float(column)
                column = f"{column:,.1f}"
            except ValueError:
                pass
        result.append(column)

    return " & ".join(result) + "\\\\"

def write_table(output_file, caption, column_count, rows, midrules = [], meta = False):
    midrules = set(midrules)
    output_file.write(table_header(caption, column_count))
    for index, row in enumerate(rows):
        if meta: row = format_metatable_row_numbers(row)
        else: row = format_row_numbers(row)
        output_file.write(row)
        if index in midrules:
            output_file.write("\\hline")
        output_file.write('\n')
    output_file.write(table_footer)

def write_image(output_file, caption, file, natwidth, natheight):
    hasher = hashlib.sha3_512()
    hasher.update(file.encode())
    hashlinkdir = os.path.join(hashdir, str(hasher.hexdigest()))
    hashlink = os.path.join(hashlinkdir, "file." + file.split(".")[-1])

    try:
        os.remove(hashlink)
    except OSError:
        pass

    pathlib.Path(hashlinkdir).mkdir(parents=True, exist_ok=True)

    try:
        os.symlink(os.path.abspath(file), hashlink)
    except FileExistsError:
        print("error: could not create symlink because it already exists. We assume that it is correct.")
    except OSError:
        print("Error: could not create symlink, but just continuing because it might have been correctly created concurrently.")
        traceback.print_exc()

    pixel_pt_factor = 0.7
    output_file.write("\\begin{figure*}\n")
    output_file.write("\\centering\n")
    output_file.write("\\includegraphics[width=\\textwidth]{." + str(hashlink) + "}\n")
    output_file.write("\\caption{" + str(caption) + "}")
    output_file.write("\\end{figure*}\n")

def convert_to_int_or_float(number):
    print("NUMBER:", number)
    try:
        return int(number)
    except ValueError:
        return float(number)

def get_values(string):
        pattern = r' ([-+]?\d+\.?\d*)? '
        numbers = re.findall(pattern, string)
        number_list = [convert_to_int_or_float(x.strip()) for x in numbers]
        return number_list


def calculate_improvement(value, comparison):
    if comparison == 0: return float('NaN')
    return value * 100 / comparison - 100

revision = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
output_file = open(output_file_name, 'w')
output_file.write(
    """
    \\documentclass[12pt,a4paper]{article}
    \\usepackage[margin=0pt]{geometry}
    \\usepackage{lmodern}
    \\usepackage{pdflscape}
    \\usepackage[T1]{fontenc}
    \\usepackage{graphicx}
    \\begin{document}
    \\begin{landscape}
    \\fontsize{6pt}{7pt}\\selectfont
    \\begin{description}
        \\item[Attention:] this file was produced automatically, and some statistics might not make sense for certain pipelines.
        \\item[Revision:] """ + str(revision) + """
    \\end{description}
    This file contains statistics about the following genome(s):
    \\begin{itemize}
    """
)

for line in name_lines:
    output_file.write("\\item " + line)
output_file.write("\\end{itemize}\n")




# Reporting graph statistics
import pandas as pd
graph_df = pd.read_csv(nodes_and_edges, sep='\t', index_col=0)
n_char = 0
with open(number_of_characters, 'r') as characters:
    n_char = characters.read()
# results = "Graph statistics & " + str(graph_df["nodes"][0]) + " & " + str(graph_df["edges"][0]) + " & " + str(graph_df["edges in cycles"][0]) + " & " + n_char  + " \\\\"
# first_line = "Parameter & nodes & edges & number of edges in all cycles & number of characters\\\\ \\hline\\\\"

first_line = "Parameter & Graph statistics\\\\ \\hline\\\\"
nodes = "Nodes in graph & " + str(graph_df["nodes"][0]) + " \\\\"
edges = "Edges in graph & " + str(graph_df["edges"][0]) + " \\\\"
edges_in_cycles = "Number of edges in all cycles & " + str(graph_df["edges in cycles"][0]) + " \\\\"
n_characters = "Number of characters in metagenome & " + n_char + " \\\\"

graph_statistics_table = [first_line, nodes, edges, edges_in_cycles, n_characters]




write_table(output_file, "Genome Graph Statistics", 1, graph_statistics_table)

write_table(output_file, "Algorithm Statistics", len(experiments), algorithm_table)

write_table(output_file, "ContigValidator", len(experiments), contig_validator_table)


# Calculating improvements
average_length = quast_table[16]
median_length = quast_table[17]
ea5max = quast_table[53]
ea10max = quast_table[54]
ea15max = quast_table[55]
ea20max = quast_table[56]
ea25max = quast_table[57]
ea30max = quast_table[58]
ea35max = quast_table[59]
ea40max = quast_table[60]
ea45max = quast_table[61]
ea50max = quast_table[62]
ea75max = quast_table[73]

integer_average_lengths = get_values(average_length)
average_length = "Improvement in average length of contigs (\%) & " + str(round(calculate_improvement(integer_average_lengths[0], integer_average_lengths[0]), 1)) + " & " + str(round(calculate_improvement(integer_average_lengths[1], integer_average_lengths[0]), 1)) + " & " + str(round(calculate_improvement(integer_average_lengths[2], integer_average_lengths[0]), 1)) + " \\\\"

integer_median_lengths = get_values(median_length)
median_length = "Improvement in median length of contigs (\%) & " + str(round(calculate_improvement(integer_median_lengths[0], integer_median_lengths[0]), 1)) + " & " + str(round(calculate_improvement(integer_median_lengths[1], integer_median_lengths[0]), 1)) + " & " + str(round(calculate_improvement(integer_median_lengths[2], integer_median_lengths[0]), 1)) + " \\\\"


integer_ea5max = get_values(ea5max)
ea5max = "Improvement in EA5max (\%) & " + str(round(calculate_improvement(integer_ea5max[0], integer_ea5max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea5max[1], integer_ea5max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea5max[2], integer_ea5max[0]), 1)) + " \\\\"
integer_ea10max = get_values(ea10max)
ea10max = "Improvement in EA10max (\%) & " + str(round(calculate_improvement(integer_ea10max[0], integer_ea10max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea10max[1], integer_ea10max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea10max[2], integer_ea10max[0]), 1)) + " \\\\"
integer_ea15max = get_values(ea15max)
ea15max = "Improvement in EA15max (\%) & " + str(round(calculate_improvement(integer_ea15max[0], integer_ea15max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea15max[1], integer_ea15max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea15max[2], integer_ea15max[0]), 1)) + " \\\\"
integer_ea20max = get_values(ea20max)
ea20max = "Improvement in EA20max (\%) & " + str(round(calculate_improvement(integer_ea20max[0], integer_ea20max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea20max[1], integer_ea20max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea20max[2], integer_ea20max[0]), 1)) + " \\\\"
integer_ea25max = get_values(ea25max)
ea25max = "Improvement in EA25max (\%) & " + str(round(calculate_improvement(integer_ea25max[0], integer_ea25max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea25max[1], integer_ea25max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea25max[2], integer_ea25max[0]), 1)) + " \\\\"
integer_ea30max = get_values(ea30max)
ea30max = "Improvement in EA30max (\%) & " + str(round(calculate_improvement(integer_ea30max[0], integer_ea30max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea30max[1], integer_ea30max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea30max[2], integer_ea30max[0]), 1)) + " \\\\"
integer_ea35max = get_values(ea35max)
ea35max = "Improvement in EA35max (\%) & " + str(round(calculate_improvement(integer_ea35max[0], integer_ea35max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea35max[1], integer_ea35max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea35max[2], integer_ea35max[0]), 1)) + " \\\\"
integer_ea40max = get_values(ea40max)
ea40max = "Improvement in EA40max (\%) & " + str(round(calculate_improvement(integer_ea40max[0], integer_ea40max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea40max[1], integer_ea40max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea40max[2], integer_ea40max[0]), 1)) + " \\\\"
integer_ea45max = get_values(ea45max)
ea45max = "Improvement in EA45max (\%) & " + str(round(calculate_improvement(integer_ea45max[0], integer_ea45max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea45max[1], integer_ea45max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea45max[2], integer_ea45max[0]), 1)) + " \\\\"
integer_ea50max = get_values(ea50max)
ea50max = "Improvement in EA50max (\%) & " + str(round(calculate_improvement(integer_ea50max[0], integer_ea50max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea50max[1], integer_ea50max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea50max[2], integer_ea50max[0]), 1)) + " \\\\"
integer_ea75max = get_values(ea75max)
ea75max = "Improvement in EA75max (\%) & " + str(round(calculate_improvement(integer_ea75max[0], integer_ea75max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea75max[1], integer_ea75max[0]), 1)) + " & " + str(round(calculate_improvement(integer_ea75max[2], integer_ea75max[0]), 1)) + " \\\\"

write_table(output_file, "QUAST: improvements of flowtigs compared to other algorithms", len(experiments), [quast_table[0]] + [average_length + "\\hline " + ea50max + ea75max + "\\hline " + ea5max + ea10max + ea15max + ea20max + ea25max + ea30max + ea35max + ea40max + ea45max], meta = True)


write_table(output_file, "QUAST: \\# of contigs", len(experiments), quast_table[0:7])
write_table(output_file, "QUAST: total length of contigs", len(experiments), [quast_table[0]] + quast_table[7:13])
write_table(output_file, "QUAST: statistics", len(experiments), [quast_table[0]] + quast_table[13:32])
write_table(output_file, "QUAST: alignment statistics", len(experiments), [quast_table[0]] + quast_table[32:], midrules = [12, 21])
write_table(output_file, "QUAST: misassembly statistics", len(experiments), quast_misassemblies_table)


# Calculating resources
import pandas as pd
resource_df = pd.read_csv(runtimes_file, sep='\t', index_col=0)
time_usage = "Runtime (s) & " + str(resource_df["time"]["unitigs"]) + " & " + str(resource_df["time"]["t. omnitigs"]) + " & " + str(resource_df["time"]["flowtigs"] + resource_df["time"]["node_to_arc"]) + " & " + str(resource_df["time"]["flowtigs"]) + "&" + str(resource_df["time"]["node_to_arc"]) + " \\\\"
memory_usage = "Memory (mb) & " + str(round(resource_df["mem"]["unitigs"]/1000, 1)) + " & " + str(round(resource_df["mem"]["t. omnitigs"]/1000, 1)) + " & " + str(round((resource_df["mem"]["flowtigs"] + resource_df["mem"]["node_to_arc"])/1000, 1)) + " & " + str(round(resource_df["mem"]["flowtigs"]/1000, 1)) + "&" + str(round(resource_df["mem"]["node_to_arc"]/1000, 1)) + " \\\\"
first_line = "Parameter & unitigs & trivial-omnitigs & flowtigs (total) & only flowtigs & only node-to-arc\\\\ \\hline\\\\"
resources_table = [first_line, time_usage, memory_usage]

write_table(output_file, "Resource usage", len(experiments) + 2, resources_table)


from os import path

if combined_eaxmax_plot_name != 'none':
    output_file.write("\\newpage")
    write_image(output_file, "EAxmax", combined_eaxmax_plot_name, 1000, 1000)

output_file.write("\\newpage")
for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
    plot_file_path = os.path.join(quast_directory, "aligned_stats", "EAxmax_plot.pdf")
    if path.isfile(plot_file_path):
        quast_png_name = plot_file_path
        write_image(output_file, "QUAST EAxmax graph for " + label, quast_png_name, 1000, 1000)
    else:
        print("Did not find '" + plot_file_path + "'")

for (label, tig_algo_file, quast_directory, contig_validator_name, resources_evaluation_name) in experiments:
    plot_file_path = os.path.join(quast_directory, "aligned_stats", "NGAx_plot.pdf")
    if path.isfile(plot_file_path):
        quast_png_name = plot_file_path
        write_image(output_file, "QUAST NGAx graph for " + label, quast_png_name, 1000, 1000)
    else:
        print("Did not find '" + plot_file_path + "'")

if bandage_png_name != 'none':
    output_file.write("\\newpage")
    write_image(output_file, "Bandage genome graph", bandage_png_name, 1000, 1000)


output_file.write(
    """
    \\end{landscape}
    \\end{document}
    """
)
