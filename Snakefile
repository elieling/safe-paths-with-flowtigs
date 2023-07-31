import itertools, sys, traceback, pathlib, os, time

###############################
###### Preprocess Config ######
###############################

print("Preprocessing config", flush = True)

# Some standard directories
# If you have a lot of data, you can change this to a larger partition.
# Keep in mind though that the $WRKDIR on turso is very fragile and you will likely shoot down the whole cluster if you do not read and understand this first:
# https://wiki.helsinki.fi/display/it4sci/Lustre+User+Guide#LustreUserGuide-4.1FileLockinginLustreisnotanadvisorytobeignored
DATADIR = "data/"
if "datadir" in config:
    DATADIR = config["datadir"]
DATADIR = os.path.abspath(DATADIR)

# If you have separate executables you can store them here.
PROGRAMDIR = "data/program"
if "programdir" in config:
    PROGRAMDIR = config["programdir"]
PROGRAMDIR = os.path.abspath(PROGRAMDIR)

# You can have a separate dir to store the plots/tables generated by your experiments.
REPORTDIR = "data/reports/"
if "reportdir" in config:
    REPORTDIR = config["reportdir"]
REPORTDIR = os.path.abspath(REPORTDIR)

print("DATADIR: {}".format(DATADIR))
print("PROGRAMDIR: {}".format(PROGRAMDIR))
print("REPORTDIR: {}".format(REPORTDIR))

# This is the maximum number of cores for your default ukko2 machine. Vorna machines have 32.
MAX_THREADS = 56
print("Setting MAX_THREADS to " + str(MAX_THREADS), flush = True)

# If you have source code in your repository, you can make snakemake compile it automatically on change by collecting all the source files.
# This is for Rust code, you can apply it to your programming language of choice.
# rust_sources = list(map(str, itertools.chain(pathlib.Path('implementation').glob('**/Cargo.toml'), pathlib.Path('implementation').glob('**/*.rs'))))

# In case you want to access the current date.
import datetime
today = datetime.date.today().isoformat()

# Chosing the number of CPUs used while running
build_cpus = 1

print("Finished config preprocessing", flush = True)

# Always use conda
workflow.use_conda = True

############################
###### File Constants ######
############################

# I like to have global constants for most involved file patterns to prevent typos.
# I keep single letter abreviations in front of the wildcards to make the resulting paths more readable.
# Helps especially when there are many arguments.
READS_SRA = os.path.join(DATADIR, "reads", "s{species}", "sra", "reads.sra")
READS_FASTA = os.path.join(DATADIR, "reads", "s{species}", "fasta", "reads.fa")
ASSEMBLY_SOURCE_READS = os.path.join(DATADIR, "assembly", "s{species}-k{k}", "reads.fa")
ASSEMBLY = os.path.join(DATADIR, "assembly", "s{species}-k{k}", "assembly.fa")
REPORT = os.path.join(REPORTDIR, "s{species}-k{k}", "report.txt")
GENOME_ALL_REFERENCES = os.path.join(DATADIR, "{file_name}.fasta")
GENOME_CONCAT_REFERENCES = os.path.join(DATADIR, "{file_name}_concat.fasta")
GENOME_CIRCULAR_REFERENCES = os.path.join(REPORTDIR, "circularization", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
# GENOME_CIRCULAR_CONCAT_REFERENCES = os.path.join(REPORTDIR, "circularization", "{file_name}_concat_k{k}ma{min_abundance}t{threads}", "report.fasta")
BUILD_FA = os.path.join(REPORTDIR, "safe_paths_unitigs", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
# BUILD_FA = os.path.join(REPORTDIR, "bcalm2", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
BUILD_LOG = os.path.join("logs", "build_{file_name}_k{k}ma{min_abundance}t{threads}", "log.log")
NODE_TO_ARC_CENTRIC_DBG_BINARY = os.path.abspath("external_software/node-to-arc-centric-dbg/target/release/node-to-arc-centric-dbg")
NODE_TO_ARC_CENTRIC_DBG = os.path.join(REPORTDIR, "node_to_arc", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.edgelist") #threads 28
SAFE_PATHS_BINARY = os.path.abspath("external_software/safe-paths/target/release/flow_decomposition")
SAFE_PATHS = os.path.join(REPORTDIR, "safe_paths_flowtigs", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta") 
# EXTERNAL_SOFTWARE_ROOTDIR = os.path.join(DATADIR, "external_software")
# QUAST_BINARY = os.path.join(EXTERNAL_SOFTWARE_ROOTDIR, "quast", "quast.py")
QUAST_BINARY = os.path.abspath("external_software/quast/quast.py")
QUAST_OUTPUT_DIR = os.path.join(REPORTDIR, "quast_{algorithm}", "{file_name}_k{k}ma{min_abundance}t{threads}")
PRACTICAL_OMNITIGS_BINARY = os.path.abspath("external_software/practical-omnitigs/implementation/target/release/cli")
PRACTICAL_OMNITIGS = os.path.join(REPORTDIR, "safe_paths_multi-safe", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
ECOLI = os.path.join(DATADIR, "ecoli.fasta")
ECOLI_CONCAT = os.path.join(DATADIR, "ecoli_concat.fasta")
SINGLE = os.path.join(DATADIR, "{file_name}.fasta")
SINGLE_CONCAT = os.path.join(DATADIR, "{file_name}_concat.fasta")
META_BASE7 = os.path.join(DATADIR, "meta_base7.fasta")
META_BASE7_CONCAT = os.path.join(DATADIR, "meta_base7_concat.fasta")
MEDIUM20 = os.path.join(DATADIR, "meta_medium20.fasta")
MEDIUM20_CONCAT = os.path.join(DATADIR, "meta_medium20_concat.fasta")
PRACTICAL_TEST_OMNITIGS = os.path.join(REPORTDIR, "safe_paths_omnitigs", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
PRACTICAL_TRIVIAL_OMNITIGS = os.path.join(REPORTDIR, "safe_paths_trivial-omnitigs", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta")
ALGORITHMS = ["unitigs", "trivial-omnitigs", "multi-safe", "flowtigs", "omnitigs"] # values for the wildcard that chooses which tigs to generate
ALGORITHM_COLUMN_NAMES = ["unitigs", "t. omnitigs", "multi-safe", "flowtigs", "omnitigs"] # column names for the different tigs
CONVERT_VALIDATION_OUTPUTS_TO_LATEX_SCRIPT = "scripts/convert_validation_outputs_to_latex.py"
CREATE_COMBINED_EAXMAX_PLOT_SCRIPT = "scripts/create_combined_eaxmax_plot.py"
REPORT_SUBDIR = os.path.join(REPORTDIR, "final_reports_{file_name}k{k}ma{min_abundance}t{threads}")
REPORT_COMBINED_EAXMAX_PLOT = os.path.join(REPORT_SUBDIR, "combined_eaxmax_plot.pdf")
REPORT_NAME_FILE = os.path.join(REPORT_SUBDIR, "name.txt")
REPORT_HASHDIR = os.path.join(REPORTDIR, "hashdir")
META_BASE7_DIR = os.path.join(DATADIR, "meta", "base7")
UNPREPROCESSED_META_BASE7 = os.path.join(META_BASE7_DIR, "GCF_0{genome}_genomic.fna")
UNPREPROCESSED_MEDIUM20 = os.path.join(DATADIR, "meta", "medium20", "GCF_{genome}_genomic.fna")
PREPROCESSED_GENOME = os.path.join(DATADIR, "metabase7", "GCF_0{genome}_genomic.fasta")
PREPROCESSED_MEDIUM = os.path.join(DATADIR, "medium20", "GCF_{genome}_genomic.fasta")
META_BASE7_UNPREPROCESSED = [os.path.join(META_BASE7_DIR, "GCF_000005845.2_ASM584v2_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_000159115.1_ASM15911v1_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_000175375.1_ASM17537v1_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_000376705.1_ASM37670v1_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_003999335.1_ASM399933v1_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_017638885.1_ASM1763888v1_genomic.fna"), os.path.join(META_BASE7_DIR, "GCF_019890915.1_ASM1989091v1_genomic.fna")]
META_BASE7_FASTA = [os.path.join(DATADIR, "metabase7", "GCF_000005845.2_ASM584v2_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_000159115.1_ASM15911v1_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_000175375.1_ASM17537v1_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_000376705.1_ASM37670v1_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_003999335.1_ASM399933v1_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_017638885.1_ASM1763888v1_genomic.fasta"), os.path.join(DATADIR, "metabase7", "GCF_019890915.1_ASM1989091v1_genomic.fasta")]
MEDIUM20_FASTA = [os.path.join(DATADIR, "medium20", "GCF_000005845.2_ASM584v2_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000009045.1_ASM904v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000013285.1_ASM1328v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000025905.1_ASM2590v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000092125.1_ASM9212v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000143845.1_ASM14384v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000159115.1_ASM15911v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000175375.1_ASM17537v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000184345.1_ASM18434v2_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000185725.2_Segn_rugo_CDC_945_V2_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000196035.1_ASM19603v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000265425.1_ASM26542v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000325705.1_ASM32570v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_000376705.1_ASM37670v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_003999335.1_ASM399933v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_008632635.1_ASM863263v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_017638885.1_ASM1763888v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_019704535.1_ASM1970453v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_019890915.1_ASM1989091v1_genomic.fasta"), os.path.join(DATADIR, "medium20", "GCF_900638215.1_56553_A01_genomic.fasta")]
META_BASE7_ABUNDANCES = os.path.join(META_BASE7_DIR, "nanosim.abundances.tsv")
MEDIUM20_ABUNDANCES = os.path.join(DATADIR, "meta", "medium20", "nanosim.abundances.tsv")
REPORT_TEX = os.path.join(REPORTDIR, "output", "{file_name}_k{k}ma{min_abundance}t{threads}", "{report_name}", "{report_file_name}.tex")



#     DATADIR = ... # wherever you have your data, e.g. /wrk-vakka/users/<your username>/flowtigs

#################################
###### Global report rules ######
#################################

# By default, snakemake executes the first rule if no target is given.
# If you don't like that, this rule makes it do nothing if no target is given.
localrules: do_nothing
rule do_nothing:
    shell:  "echo 'No target specified'"

# Here can be rules that specify a set of reports you would like to create, like e.g. all executions for all combinations of different k values and species.
localrules: report_all
rule report_all:
    input:  reports = expand(REPORT, species = ["ecoli", "scerevisiae"], k = [23, 31, 39]),

###############################
###### Report Generation ######
###############################

# Here can be rules to generate reports from your experiments, like e.g. plots or tables.
localrules: create_single_report
rule create_single_report:
    input:  assembly = ASSEMBLY,
    output: report = REPORT,
    conda:  "config/conda-biopython-env.yml"
    script: "scripts/create_single_report.py"


class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

def safe_format(str, **kwargs):
    return str.format_map(SafeDict(kwargs))

def safe_expand(str, **kwargs):
    items = []
    for key, values in kwargs.items():
        if type(values) is str or type(values) is not list:
            values = [values]
        items.append([(key, value) for value in values])

    for combination in itertools.product(*items):
        yield safe_format(str, **dict(combination))

def wildcard_format(str, wildcards):
    return str.format(**dict(wildcards.items()))


def get_single_report_script_column_arguments_from_wildcards(wildcards):
    try:
        result = ""
        once = True
        for algorithm in ALGORITHMS:
            if once:
                once = False
            else:
                result += " "

            quast_output_dir = safe_format(QUAST_OUTPUT_DIR, algorithm = algorithm).format(**wildcards)
            result += f"'{algorithm}' '' '{quast_output_dir}' '' ''"
        return result
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")



rule create_single_report_tex:
    input:  quasts = [safe_format(QUAST_OUTPUT_DIR, algorithm = algorithm) for algorithm in ALGORITHMS],
            combined_eaxmax_plot = REPORT_COMBINED_EAXMAX_PLOT,
            script = CONVERT_VALIDATION_OUTPUTS_TO_LATEX_SCRIPT,
    output: report = REPORT_TEX,
    log:    log = "logs/create_single_report/{file_name}_k{k}ma{min_abundance}t{threads}r{report_name}rf{report_file_name}/log.log",
    params: genome_name = lambda wildcards: ", ".join(wildcards.file_name), 
            script_column_arguments = get_single_report_script_column_arguments_from_wildcards,
            name_file = REPORT_NAME_FILE,
            hashdir = REPORT_HASHDIR,
    wildcard_constraints:
            report_name = "[^/]+",
    conda: "config/conda-latex-gen-env.yml" 
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        mkdir -p '{params.hashdir}'
        echo '{wildcards.report_name} {params.genome_name} {wildcards.report_file_name}' > '{params.name_file}'
        python3 '{input.script}' '{params.hashdir}' '{params.name_file}' 'none' 'none' '{input.combined_eaxmax_plot}' '{output}' {params.script_column_arguments}
        """


rule create_combined_eaxmax_graph:
    input:  quast_csvs = [os.path.join(safe_format(QUAST_OUTPUT_DIR, algorithm = algorithm), "aligned_stats", "EAxmax_plot.csv") for algorithm in ALGORITHMS],
            script = CREATE_COMBINED_EAXMAX_PLOT_SCRIPT,
    output: REPORT_COMBINED_EAXMAX_PLOT,
    params: input_quast_csvs = lambda wildcards, input: "' '".join([shortname + "' '" + quast for shortname, quast in zip(ALGORITHM_COLUMN_NAMES, input.quast_csvs)])
    conda:  "config/conda-seaborn-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        mkdir -p "$(dirname '{output}')"
        python3 '{input.script}' '{params.input_quast_csvs}' '{output}'
        """

    

########################
###### Assembly ########
########################


# Rule to get sequences for metagenome data.
localrules: metagenome_to_single_file
rule metagenome_to_single_file:
    input:  abundances = META_BASE7_ABUNDANCES,
            references = META_BASE7_FASTA,
    output: report = META_BASE7,
    log:    log = "logs/metagenome_to_single_file/log.log",
    conda:  "config/conda-seaborn-env.yml"
    script: "scripts/metagenome_to_single_file.py"



# Rule to concatenate metagenome data.
localrules: metagenome_concatenate
rule metagenome_concatenate:
    input:  references = META_BASE7_FASTA,
    output: report = META_BASE7_CONCAT,
    log:    log = "logs/metagenome_concatenate/log.log",
    conda:  "config/conda-seaborn-env.yml"
    script: "scripts/metagenome_concatenate.py"


# Rule to get sequences for metagenome data.
localrules: medium20_to_single_file
rule medium20_to_single_file:
    input:  abundances = MEDIUM20_ABUNDANCES,
            references = MEDIUM20_FASTA,
    output: report = MEDIUM20,
    log:    log = "logs/medium20_to_single_file/log.log",
    conda:  "config/conda-seaborn-env.yml"
    script: "scripts/metagenome_to_single_file.py"



# Rule to concatenate metagenome data.
localrules: medium20_concatenate
rule medium20_concatenate:
    input:  references = MEDIUM20_FASTA,
    output: report = MEDIUM20_CONCAT,
    log:    log = "logs/medium20_concatenate/log.log",
    conda:  "config/conda-seaborn-env.yml"
    script: "scripts/metagenome_concatenate.py"



# Rule to remove all the non-[A, C, G, T] characters from a fasta file. Should be ran first.
rule preprocessing_single_genome:
    input:  assembly = UNPREPROCESSED_META_BASE7,
    output: report = PREPROCESSED_GENOME,
    log:    log = "logs/preprocessing_single_genome/{genome}/log.log",
    conda:  "config/conda-biopython-env.yml"
    script: "scripts/preprocessing.py"


# Rule to remove all the non-[A, C, G, T] characters from a fasta file. Should be ran first.
rule preprocessing_single_medium20_genome:
    input:  assembly = UNPREPROCESSED_MEDIUM20,
    output: report = PREPROCESSED_MEDIUM,
    log:    log = "logs/preprocessing_single_genome/{genome}/log.log",
    conda:  "config/conda-biopython-env.yml"
    script: "scripts/preprocessing.py"


# Rule to circularize non-circular sequences of the genome or metagenome.
# input: data of the non-circular sequences of the genome or metagenome in fasta format.
# wildcards: k=size of kmers, file_name=name of the file with the data in the data folder,
#   min_abundance=minimum abundance, threds=number of cpu cores used.
# output: data of the circular sequences of the genome or metagenome in fasta format.
localrules: circularization
rule circularization:
    input:  assembly = GENOME_ALL_REFERENCES,
    output: report = GENOME_CIRCULAR_REFERENCES,
    log:    log = "logs/circularization/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    conda:  "config/conda-biopython-env.yml"
    script: "scripts/circularization.py"


# Rule to create a node node-cenric weighted De Bruijn graph using the bcalm2 tool.
# input: data of the circular sequences of the genome or metagenome in fasta format.
# wildcards: k=size of kmers, file_name=name of the file with the data in the data folder,
#   min_abundance=minimum abundance, threds=number of cpu cores used.
# output: Node-centric weighted De Bruijn graph with k-mers of size k in fasta format.
rule bcalm2_build:
    input:  references = GENOME_CIRCULAR_REFERENCES,
    output: tigs = BUILD_FA,
    log:    log = "logs/bcalm2/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    params: references = lambda wildcards, input: "'" + "' '".join(input.references) + "'",
    conda:  "config/conda-bcalm2-env.yml",
    threads: MAX_THREADS,
    shadow: "minimal"
    resources:
            mem_mb = 10_000, # probably much more than needed
            time_min = 60,
            cpus = build_cpus,
            queue = 'aurinko,bigmem,short,medium', # I had some more complex expression here, the queues fitting to the time are on https://wiki.helsinki.fi/display/it4sci/HPC+Environment+User+Guide#HPCEnvironmentUserGuide-4.8.4Partitions-Ukko
    shell:  """
        rm -f '{log.log}'
        ${{CONDA_PREFIX}}/bin/time -v bcalm -nb-cores {threads} -kmer-size {wildcards.k} -in '{input.references}' -out '{output.tigs}' -abundance-min {wildcards.min_abundance} 2>&1 | tee -a '{log.log}'
        mv '{output.tigs}.unitigs.fa' '{output.tigs}'
        """


# Rule to set up the external software for the node_to_arc_centric_dbg rule.
rule build_node_to_arc_centric_dbg:
    input:  "external_software/node-to-arc-centric-dbg/Cargo.toml",
    output: NODE_TO_ARC_CENTRIC_DBG_BINARY,
    conda:  "config/conda-rust-env.yml",
    threads: MAX_THREADS,
    resources:
            mem_mb = 10000,
            time_min = 60,
            cpus = MAX_THREADS,
            queue = "aurinko,bigmem,short,medium",
    shell:  """
        cd external_software/node-to-arc-centric-dbg
        cargo build --release -j {threads} --offline
    """


# Rule to transform the node-centric weighted De Bruijn graph given by the bcalm2_build
#   rule to an arc-centric weighted De Bruijn graph.
# input: output of bcalm2_build and output of build_node_to_arc_centric_dbg.
# wildcards: k=size of kmers, file_name=name of the file with the data in the data folder,
#   min_abundance=minimum abundance, threads=number of cpu cores used.
# output: Arc-centric weighted De Bruijn graph with k-mers of size k in fasta format.
rule node_to_arc_centric_dbg:
    input:  node_centric_dbg = BUILD_FA,
            binary = NODE_TO_ARC_CENTRIC_DBG_BINARY,
    log:    log = "logs/node_to_arc/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    output: arc_centric_dbg = NODE_TO_ARC_CENTRIC_DBG,   # .edgelist
    conda:  "config/conda-time-env.yml",
    resources:
            time_min = 60, # likely too much
            mem_mb = 10_000, # likely too much
            queue = "short,medium,bigmem,aurinko",
    shell:  """
        rm -f '{log.log}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -k {wildcards.k} --input '{input.node_centric_dbg}' --output '{output.arc_centric_dbg}'
    """




# Rule to set up the external software for the safe-paths rule.
rule build_safe_paths:
    input:  "external_software/safe-paths/Cargo.toml",
    output: SAFE_PATHS_BINARY,
    conda:  "config/conda-rust-env.yml",
    threads: MAX_THREADS,
    resources:
            mem_mb = 10000,
            time_min = 60,
            cpus = MAX_THREADS,
            queue = "aurinko,bigmem,short,medium",
    shell:  """
        cd external_software/safe-paths
        cargo build --release -j {threads} --offline
    """


# Rule to get the safe paths from the edge-centric weighted De Bruijn graph given by
#   the node_to_arc_centric_dbg rule.
# input: output of node_to_arc_centric_dbg and output of build_safe_paths.
# wildcards: k=size of kmers, file_name=name of the file with the data in the data folder,
#   min_abundance=minimum abundance, thredas=number of cpu cores used.
# output: List of the safe paths with, for each path, a line with ">Path_<path-index>"
#   and a second line with the sequence.
rule safe_paths:
    input:  arc_centric_dbg = NODE_TO_ARC_CENTRIC_DBG,
            binary = SAFE_PATHS_BINARY,
    log:    log = "logs/safe_paths/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    output: safe_paths = SAFE_PATHS,
    conda:  "config/conda-time-env.yml",
    resources:
            time_min = 60, # likely too much
            mem_mb = 10_000, # likely too much
            queue = "short,medium,bigmem,aurinko",
    shell:  """
        rm -f '{log.log}'
        ${{CONDA_PREFIX}}/bin/time -v '{input.binary}' -k {wildcards.k} --input '{input.arc_centric_dbg}' --output '{output.safe_paths}'
    """




# Rule to make a report of the results.
# input: 
#   contigs: List of safe paths in fasta format.
#   script: File to execute the Rust program.
#   error_report: Template file for when there is an error creating the report.tex file.
#   error_missassemblies_report: Template file for when there is an error creating the missassemblies_report.tex file.
# wildcards: 
#   algorithm: method for finding safe paths. Either "uni" for unitigs, "omni" for omnitigs or "flow" for flowtigs.
#   k: size of kmers
#   file_name: name of the file with the data in the data folder
#   min_abundance: minimum abundance
#   threads: number of cpu cores used.
# output: report in .tex, .tsv, .txt, and .pdf formats
rule run_quast:
    input:  contigs = os.path.join(REPORTDIR, "safe_paths_{algorithm}", "{file_name}_k{k}ma{min_abundance}t{threads}", "report.fasta"),
            references = [GENOME_CONCAT_REFERENCES], # list of references
            script = QUAST_BINARY,
            error_report = "config/quast_error_report.tex",
            error_misassemblies_report = "config/quast_error_misassemblies_report.tex"
    output: directory = directory(QUAST_OUTPUT_DIR),
            eaxmax_csv = os.path.join(QUAST_OUTPUT_DIR, "aligned_stats/EAxmax_plot.csv"),
    log:    minimap = os.path.join(QUAST_OUTPUT_DIR, "contigs_reports", "contigs_report_contigs.stderr")
    params: references = lambda wildcards, input: "-r '" + "' -r '".join(input.references) + "'",
    conda: "config/conda-quast-env.yml"
    threads: 14,
    resources: mem_mb = 10_000, # likely to much for our genomes
               cpus = 14,
               time_min = 60,
               queue = "short,medium,bigmem,aurinko", # I had some more complex expression here, the queues fitting to the time are on https://wiki.helsinki.fi/display/it4sci/HPC+Environment+User+Guide#HPCEnvironmentUserGuide-4.8.4Partitions-Ukko
    shell:  """
        set +e 
        ${{CONDA_PREFIX}}/bin/time -v {input.script} -t {threads} --no-html -o '{output.directory}' {params.references} '{input.contigs}'
        set -e
        
        if [ $? -ne 0 ]; then
            rm -rf '{output.directory}'
        fi


        mkdir -p '{output.directory}'
        if [ ! -f '{output.directory}/report.tex' ]; then
            echo "report.tex is missing, using error template" 
            cp '{input.error_report}' '{output.directory}/report.tex'
        fi

        mkdir -p '{output.directory}/contigs_reports'
        if [ ! -f '{output.directory}/contigs_reports/misassemblies_report.tex' ]; then
            echo "misassemblies_report.tex is missing, using error template"
            cp '{input.error_misassemblies_report}' '{output.directory}/contigs_reports/misassemblies_report.tex'
        fi

        mkdir -p '{output.directory}/aligned_stats'
        if [ ! -f '{output.directory}/aligned_stats/EAxmax_plot.csv' ]; then
            echo "EAxmax_plot.csv is missing, using empty file"
            touch '{output.directory}/aligned_stats/EAxmax_plot.csv'
        fi
    """ #  --large




rule latex:
    input: "{subpath}{file_suffix}.tex"
    output: "{subpath}{file_suffix}.pdf"
    wildcard_constraints:
        file_suffix = "((report)|(all)|(reduced))",
    conda: "config/conda-latex-env.yml"
    threads: 1
    resources:
            queue = "short,medium,bigmem,aurinko",
    shell: """
        tectonic '{input}'
        """




# Rule to set up the external software for the practical omnitigs rule.
rule build_practical_omnitigs:
    input:  "external_software/practical-omnitigs/implementation/Cargo.toml",
    output: PRACTICAL_OMNITIGS_BINARY,
    conda:  "config/conda-practical-omnitigs-env.yml",
    threads: MAX_THREADS,
    resources:
            mem_mb = 10000,
            time_min = 60,
            cpus = MAX_THREADS,
            queue = "aurinko,bigmem,short,medium",
    shell:  """
        cd external_software/practical-omnitigs/implementation
        cargo build --release -j {threads} 
    """



# Rule to get the safe paths from the node-centric weighted De Bruijn graph given by
#   the bcalm2 rule.
# input: output of bcalm2 and output of build_practical_omnitigs.
# wildcards: k=size of kmers, file_name=name of the file with the data in the data folder,
#   min_abundance=minimum abundance, thredas=number of cpu cores used.
# output: List of the safe paths with, for each path, a line with "><path-index>"
#   and a second line with the sequence.
rule practical_omitigs:
    input:  practical_omnitigs = BUILD_FA,
            binary = PRACTICAL_OMNITIGS_BINARY,
    log:    log = "logs/practical_omnitigs/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    output: practical_omnitigs = PRACTICAL_OMNITIGS,  
    conda:  "config/conda-rust-env.yml",
    resources:
            time_min = 60, # likely too much
            mem_mb = 10_000, # likely too much
            queue = "short,medium,bigmem,aurinko",
    shell:  """
        rm -f '{log.log}'
        '{input.binary}' compute-multi-safe-walks --file-format bcalm2 --input '{input.practical_omnitigs}' --output '{output.practical_omnitigs}' -k {wildcards.k}
    """


# Rule for testing the practical omnitigs rule
rule practical_test_omitigs:
    input:  practical_omnitigs = BUILD_FA,
            binary = PRACTICAL_OMNITIGS_BINARY,
    log:    log = "logs/practical_test_omnitigs/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    output: practical_omnitigs = PRACTICAL_TEST_OMNITIGS,  
    conda:  "config/conda-rust-env.yml",
    resources:
            time_min = 60, # likely too much
            mem_mb = 10_000, # likely too much
            queue = "short,medium,bigmem,aurinko",
    shell:  """
        rm -f '{log.log}'
        '{input.binary}' compute-omnitigs --file-format bcalm2 --input '{input.practical_omnitigs}' --output '{output.practical_omnitigs}' -k {wildcards.k}
    """


# Rule for testing the practical omnitigs rule
rule practical_trivial_omitigs:
    input:  practical_omnitigs = BUILD_FA,
            binary = PRACTICAL_OMNITIGS_BINARY,
    log:    log = "logs/practical_trivial_omnitigs/{file_name}_k{k}ma{min_abundance}t{threads}/log.log",
    output: practical_omnitigs = PRACTICAL_TRIVIAL_OMNITIGS,  
    conda:  "config/conda-rust-env.yml",
    resources:
            time_min = 60, # likely too much
            mem_mb = 10_000, # likely too much
            queue = "short,medium,bigmem,aurinko",
    shell:  """
        rm -f '{log.log}'
        '{input.binary}' compute-trivial-omnitigs --non-scc --file-format bcalm2 --input '{input.practical_omnitigs}' --output '{output.practical_omnitigs}' -k {wildcards.k}
    """



#######################
###### Downloads ######
#######################

# Here can be rules to download your input files.
# These must be localrules, as the compute nodes of the federated cluster do not have a connection to the internet.

def get_reads_url(wildcards):
    # Snakemake just swallows whatever errors happen in functions.
    # The only way around that is to handle and print the error ourselves.
    try:
        species = wildcards.species

        if species == "ecoli":
            return "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR12793243/SRR12793243.1"
        elif species == "scerevisiae":
            return "https://sra-download.ncbi.nlm.nih.gov/traces/sra78/SRR/014398/SRR14744387"
        else:
            raise Exception("Unknown species: {}".format(species))
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")


# Rule to download the E.coli bacteria data
localrules: download_ecoli
rule download_ecoli:
    output: ECOLI
   # conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p data
        cd data

        rm -rf ecoli.fasta
        wget -O ecoli.fasta.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
        gunzip ecoli.fasta.gz
    """

# Rule needed for run_quast. Just copies the file with a new name.
rule ecoli_concat:
    input: file = ECOLI
    output: file = ECOLI_CONCAT
    shell: """
        cp {input.file} {output.file}
    """

# Rule to create the the concatenated version of a single file. Just copies the file with a new name.
rule single_concat:
    input: file = SINGLE
    output: file = SINGLE_CONCAT
    shell: """
        cp {input.file} {output.file}
    """


# Rule to download the external software used for turning the node-centric De Bruijn graph
#   to an arc-centric one.
localrules: download_node_to_arc_centric_dbg
rule download_node_to_arc_centric_dbg:
    output: "external_software/node-to-arc-centric-dbg/Cargo.toml"
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p external_software
        cd external_software

        rm -rf node-to-arc-centric-dbg
        git clone https://github.com/sebschmi/node-to-arc-centric-dbg.git
        cd node-to-arc-centric-dbg
        git checkout 9ae1e61bd089be227342d03762031381cb989888 
        cargo fetch
    """ # d0afdf0532fcaa658f57182fe4e5d26224136285



# Rule to download the external software used for calculatings safe paths from an arc-centric De Bruijn graph.
localrules: download_safe_paths
rule download_safe_paths:
    output: "external_software/safe-paths/Cargo.toml"
    conda:  "config/conda-rust-env.yml"
    threads: 1
    shell:  """
        mkdir -p external_software
        cd external_software

        rm -rf safe-paths
        git clone https://github.com/elieling/safe-paths.git
        cd safe-paths
        git checkout 19304e99283e336ebf15d5b206203da00b053f82

        cargo fetch
    """ # db71ee631c4a2f5679ea1ecb06052e233a987e26 


localrules: install_quast
rule install_quast:
    output: script = QUAST_BINARY,
    params: external_software_dir = "external_software", # EXTERNAL_SOFTWARE_ROOTDIR,
    threads: 1
    shell: """
        mkdir -p {params.external_software_dir}
        cd {params.external_software_dir}

        rm -rf quast
        git clone https://github.com/elieling/quast.git 
        cd quast
        git checkout af809a144d007f8a31cfffee994b33362a3e379e
    """ # 194f0026b54849b78fc145e7b9152cecb2ebd00d f0c416b5892fc7b7caf4ba1b82bdca85761fc952 0da05f6ef5c6503dc95ab2f8b7a4a3da21cc5172



localrules: download_practical_omnitigs
rule download_practical_omnitigs:
    output: "external_software/practical-omnitigs/implementation/Cargo.toml"
    conda:  "config/conda-practical-omnitigs-env.yml"
    threads: 1
    shell:  """
        mkdir -p external_software
        cd external_software

        rm -rf practical-omnitigs
        git clone https://github.com/algbio/practical-omnitigs.git
        cd practical-omnitigs
        git checkout 690a836ebc158b17d68593fdcc59c1a05cc07b20 
        cd implementation
        cargo fetch
    """ # e1640a6220dedb85ad90f83d9b141f4fd257a6e5 bb1de69873c6b48f183e51bca2f48d2a057b8b64
        # git checkout sebschmi/unspecified 






localrules: download_sra_file
rule download_sra_file:
    output: file = READS_SRA,
    params: url = get_reads_url,
    shell:  """
        wget --progress=dot:mega -O '{output.file}' '{params.url}'
    """

# If no resources are specified, the rule will use the default resources specified in config/turso/config.yml
rule convert_sra_download:
    input:  file = READS_SRA,
    output: file = READS_FASTA,
    conda:  "config/conda-convert-reads-env.yml"
    shell:  "fastq-dump --stdout --fasta default '{input.file}' > '{output.file}'"

# I link the reads to each assembly directory separately, such that bcalm creates its auxiliary files in that directory instead of the download directory.
localrules: download_reads
rule download_reads:
    input:  file = READS_FASTA,
    output: file = ASSEMBLY_SOURCE_READS,
    shell:  "ln -sr -T '{input.file}' '{output.file}'"

##############################
###### Download results ######
##############################

# Here is a template for a rule to download report files from turso to your local machine.
# You can use something similar, just make sure to update the paths and the include specification.
#localrules: sync_turso_results
#rule sync_turso_results:
#    conda: "config/conda-rsync-env.yml"
#    shell: """
#        mkdir -p data/reports
#        rsync --verbose --recursive --no-relative --include="*/" --include="report.pdf" --include="aggregated-report.pdf" --exclude="*" turso:'/proj/sebschmi/git/practical-omnitigs/data/reports/' data/reports
#        """

