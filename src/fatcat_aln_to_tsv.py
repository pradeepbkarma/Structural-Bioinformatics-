import pandas as pd

# read the fatcat alignment file 
aln_file = "../Data/

# make list of all the columns of interests
query=[]
qlen=[]
target = []
tlen = []
opt_rmsd = []
score = []
length = [] # alignment length 
p_val = []
ident = []
similarity = []

# read each line and append appropriate information to each list 
with open(aln_file, 'r')as infile:
    for line in infile:
        line = line.strip()
        if line.startswith("Align"):
            parts = line.split()
            query.append(parts[1].split('.')[0]) 
            qlen.append(int(parts[2])) # qury length
            target.append(parts[4].split('.')[0])
            tlen.append(int(parts[5])) # Target length 
        
        elif line.startswith("Twists"):
            parts = line.split()
            rmsd_index = parts.index("opt-rmsd") + 1
            score_index = parts.index("Score") + 1
            length_index = parts.index("align-len") + 1
            opt_rmsd.append(float(parts[rmsd_index]))
            score.append(float(parts[score_index]))
            length.append(float(parts[length_index]))

        elif line.startswith("P-value"):
            parts = line.split()
            p_val_index = parts.index("P-value") + 1
            ident_index = parts.index("Identity") + 1
            sim_index = parts.index("Similarity") + 1
            p_val.append(float(parts[p_val_index]))
            ident.append(float(parts[ident_index].replace("%", "")))
            similarity.append(float(parts[sim_index].replace("%", "")))

# create a Dataframe from the extracted data

df = pd.DataFrame({
    "query": query,
    "target": target,
    "qlen": qlen,
    "tlen": tlen,
    "optimal_rmsd": opt_rmsd,
    "Score": score,
    "aln length": length,
    "p-val": p_val,
    "Identity (%)": ident,
    "similarity (%)": similarity
})

output_file = "../Data/test_allpairs_aln.tsv"
df.to_csv(output_file, sep="\t", index=False)
