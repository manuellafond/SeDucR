import re

file = open("rt_segdup.txt")
input = file.readline()
file.close()

# parse input for segrec
srSpecies = open("rt-sp.txt", "w")
srSpecies.write(re.sub(r"\w+", r"'\g<0>'", input[4 : input.find('"', 5)] + ";"))
srSpecies.close()

geneNum = 1
srGenes = open("rt-genes.txt", "w")
genesList = re.sub(" -G ", "\n", input[input.find('"', 5) + 1 :]).split("\n")[1:-1]
for g in genesList:
    g = g[1 : g.find('"', 1)]
    g = re.sub(r"\w+", r"\g<0>_" + f"{geneNum}", g) + ";"
    srGenes.write(g + "\n")
    geneNum += 1
srGenes.close()
