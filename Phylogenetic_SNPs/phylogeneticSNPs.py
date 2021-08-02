import pandas as pd
import numpy as np

matrix = pd.read_csv("merged_bovisorygiscaprae_840_filtered_filGATK_mac7_matrix",
                     sep = "\t")
# 1 corresponds to 0/0 which is the reference, 0 corresponds to missing and 3 to 1/1 which is the alternative
# thus replace it in that way
matrix = matrix.replace(0, np.nan)
matrix = matrix.replace(1, 0)
matrix = matrix.replace(3, 1)
info = pd.read_csv("metadata_bovisorygiscaprae_840.txt", 
                  sep = "\t")
info["Lineage"] = np.where(info.ClonalComplex == "orygis", "orygis", np.where(info.ClonalComplex == "caprae", "caprae", np.where(info.ClonalComplex == "L6", "L6","bovis")))
matrix = matrix.merge(info[["G_NUMBER", "ClonalComplex", "Lineage"]], left_index = True, right_on = "G_NUMBER")
# for each lineage and sublineage interested in, find the positions that are mutated in all the strains
snp_dic_lin = {}
matrix = matrix[matrix.Lineage != "L6"] # exclude the outgroup
for lineage in matrix.Lineage.unique():
    subset = matrix[matrix.Lineage == lineage]
    pos = list(subset.loc[:, (subset==1).all()].columns)
    snp_dic_lin[lineage] = pos
snp_dic_sublin = {}
for lineage in [ 'unk2', 'af2', 'af1', 'PZA_sus', 'unk3', 'eu1',
       'unk4', 'eu2', 'unk7', 'unk5', 'unk6', 'unk9', 'La1.2 BCG']:
    subset = matrix[matrix.ClonalComplex.isin([lineage])]
    pos = list(subset.loc[:, (subset==1).all()].columns)
    snp_dic_sublin[lineage] = pos
# now check, which of the positions are uniquely mutated in that lineage/sublineage
unique_snps = {}
for lineage in matrix.Lineage.unique():
    subdic = [value for key, value in snp_dic_lin.items() if key not in [lineage]]
    pos_rest = [item for sublist in subdic for item in sublist]
    pos_unique = np.setdiff1d(snp_dic_lin[lineage], pos_rest)
    unique_snps[lineage] = pos_unique
# also add all the combinations:
for combination in [["af2", "unk2"], ["unk4", "unk5", "eu2"], ["eu1", "unk6", "unk7"], ["unk2", "La1.2 BCG"]]:
    
    subdic = [value for key, value in snp_dic_sublin.items() if key not in combination]
    pos_rest = [item for sublist in subdic for item in sublist]
    combi_dic = [value for key, value in snp_dic_sublin.items() if key in combination]
    all_combi_pos = []
    for i in range(len(combination)):
        all_combi_pos.append(snp_dic_sublin[combination[i]])
    pos_combi = list(set.intersection(*map(set, all_combi_pos)))
    pos_unique = np.setdiff1d(pos_combi, pos_rest)
    name = "_".join(combination)
    unique_snps[name] = pos_unique
# add the combi bovis caprae
for combination in [["bovis", "caprae"]]:
    
    subdic = [value for key, value in snp_dic_lin.items() if key not in combination]
    pos_rest = [item for sublist in subdic for item in sublist]
    combi_dic = [value for key, value in snp_dic_lin.items() if key in combination]
    all_combi_pos = []
    for i in range(len(combination)):
        all_combi_pos.append(snp_dic_lin[combination[i]])
    pos_combi = list(set.intersection(*map(set, all_combi_pos)))
    pos_unique = np.setdiff1d(pos_combi, pos_rest)
    name = "_".join(combination)
    unique_snps[name] = pos_unique
# add the sublineages you are interested in:
snp_dic_sublin.update(snp_dic_lin) # also add the lineage-specific SNPs
for lineage in ['orygis', 'unk2', 'af2', 'af1', 'PZA_sus', 'unk3', 'eu1',
       'caprae', 'unk4', 'eu2', 'unk7', 'unk5', 'unk6', 'unk9', 'La1.2 BCG']:
    subdic = [value for key, value in snp_dic_sublin.items() if key not in [lineage]]
    pos_rest = [item for sublist in subdic for item in sublist]
    pos_unique = np.setdiff1d(snp_dic_sublin[lineage], pos_rest)
    unique_snps[lineage] = pos_unique
pos_dic = {}
i = 0
for key in list(unique_snps.keys()):
    for item in list(unique_snps[key]):
        pos_dic[i] = {"Position_ref" : item.split(":")[1].split("_")[0],
                     "ancestral" : item.split(":")[1].split("_")[1].split("/")[0],
                     "derived" : item.split(":")[1].split("_")[1].split("/")[1],
                     "PhylogeneticSNP": key}
        i += 1


# now check in what gene the positions are located and what the gene based position is
annot = pd.read_csv("AnnotationH37Rv.ptt",
                   sep = " ")
annot[["Start", "Stop"]] = annot[["Start", "Stop"]].apply(pd.to_numeric)

for i in range(len(pos_dic)):
    pos = int(pos_dic[i]["Position_ref"])
    sub = annot[(annot.Stop >= pos) & (annot.Start <= pos)]
    if len(sub) == 1:
        gene = sub.ID.item()
        strand = sub.Strand.item()
        if strand == "+":
            gene_start = sub.Start.item()
            gene_end = sub.Stop.item()
            pos_gene = int(pos) - int(gene_start)+1
            rest = pos_gene%3
            if rest == 1:
                position_codon = "first"
            elif rest == 0:
                position_codon = "third"
            elif rest == 2:
                position_codon = "second"
               
        elif strand == "-": 
            gene_end = sub.Start.item()
            gene_start = sub.Stop.item()
            pos_gene = int(gene_start) - int(pos) + 1
            rest = pos_gene%3
            if rest == 1:
                position_codon = "first"
            elif rest == 0:
                position_codon = "third"

            elif rest == 2:
                position_codon = "second"

        else:
            gene_start = sub.Start.item()
            gene_end = sub.Stop.item()
            pos_gene = int(pos) - int(gene_start)
            mutation = np.nan
            codon_new = np.nan
            codon = np.nan
            position_codon = np.nan
        pos_dic[i]["Start"] = gene_start
        pos_dic[i]["End"] = gene_end
        pos_dic[i]["Strand"] = strand
        pos_dic[i]["Gene"] = gene
        pos_dic[i]["Position_gene"] = pos_gene
        pos_dic[i]["Position_codon"] = position_codon
    else: 
        pos_dic[i]["Start"] = np.nan
        pos_dic[i]["End"] = np.nan
        pos_dic[i]["Strand"] = np.nan
        pos_dic[i]["Gene"] = np.nan
        pos_dic[i]["Position_gene"] = np.nan
        pos_dic[i]["Position_codon"] = np.nan

phylogenetic_snps = pd.DataFrame(pos_dic).T

# now add the annotations to the file containing all the phylogenetic SNPs

annot = pd.read_csv("annotations_bovisorygiscaprae_804.txt",
                   sep = " ", header = None)
dic = pd.DataFrame.to_dict(phylogenetic_snps.T)
for i in range(len(dic)):
    pos = dic[i]["Position_ref"]
    snp_eff_annot = annot[annot[1] == pos][4].item()
    mut = snp_eff_annot.split("|")[1]
    if mut == "missense_variant":
        mut = "nonsynonymous"
    elif mut == "synonymous_variant":
        mut = "synonymous"
    dic[i]["snpEff_mutation"] = mut
phylogenetic_snps = pd.DataFrame(dic).T
# check for all the lifestock-associated phylogenetic SNPs whether they are also found variable in at least one of the Treemmer
# dataset. If they are, then remove them.

pos_treemmer = pd.read_csv("merged_treemmed_mtbc_4742_filtered_filGATK_variablepos_refaltposition.txt",
                          sep = " ", header = None)
pos_treemmer.head()
homoplasic = []
for item in phylogenetic_snps.Position_ref.unique():
    if item in pos_treemmer[0].unique():
        homoplasic.append(item)
print(len(phylogenetic_snps))
red = phylogenetic_snps[~phylogenetic_snps.Position_ref.isin(homoplasic)]
print(len(red))


kvarq = pd.read_csv("subset5_phylogenetic_SNPs_bovisorygiscaprae_840_nohomoplasicpos.txt",
                     sep = "\t", header = None)
red["KvarQ_informative"] = np.where(red.Position_ref.isin(kvarq[1].unique()), True, False)

# exclude intergenic positions
red = red[red.Strand.isin(["+", "-"])]

# convert to the new naming_sceme
translation = {"orygis":"La3","bovis":"La1","caprae":"La2","unk4_unk5_eu2":"La1.7",
               "eu1_unk6_unk7":"La1.8","unk2":"La1.2","af2":"La1.3",
               "af1":"La1.6","PZA_sus":"La1.1","unk3":"La1.4","eu1":"La1.8.1",
               "unk4":"unk4","eu2":"La1.7.1","unk7":"La1.8.2","unk5":"unk5",
               "unk6":"unk6","unk9":"La1.5", 'af2_unk2': "La1.3_La1.2",
              'bovis_caprae': "La1_La2", 'unk2_La1.2 BCG':'La1.2_La1.2 BCG',
              "La1.2 BCG": "La1.2 BCG"}
red = red.reset_index(drop = True)
dic = pd.DataFrame.to_dict(red.T)
for i in range(len(dic)):
    old = dic[i]["PhylogeneticSNP"]
    new = translation[old]
    dic[i]["PhylogeneticSNP"] = new
snps = pd.DataFrame(dic).T
snps[["Start", "End", "Position_gene"]] = snps[["Start", "End", "Position_gene"]].astype(int)
snps.to_csv("phylogenetic_SNPs_bovisorygiscaprae_840_nohomoplasicpos_nointergenic_kvarqmarked.txt",
                     sep = "\t", index = False)