import sys
import pandas as pd

folder_name = sys.argv[1]
abundances_file = sys.argv[2]

abundances_df = pd.read_csv(abundances_file, sep='\t')
names = abundances_df["Size"]

result = "["
first = True
for name in names:
    if name.endswith(".fna") or name.endswith(".fasta"):
        if not first: result = result + ", "
        result = result + "os.path.join(DATADIR, '" + folder_name + "')"
        first = False

result = result + "]"

print(result)