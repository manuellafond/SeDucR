import re

guigoFile = open("guigo_all.txt")
guigoInput = guigoFile.readline()
guigoFile.close()

# parse input for segrec
srSpecies = open("guigo-sp.txt", "w")
srSpecies.write(
    re.sub(r"\w+", r"'\g<0>'", guigoInput[4 : guigoInput.find('"', 5)] + ";")
)
srSpecies.close()

geneNum = 1
srGenes = open("guigo-genes.txt", "w")
genesList = re.sub(" -G ", "\n", guigoInput[guigoInput.find('"', 5) + 1 :]).split("\n")[
    1:-1
]
for g in genesList:
    g = g[1 : g.find('"', 1)]
    g = re.sub(r"\w+", r"\g<0>_" + f"{geneNum}", g) + ";"
    srGenes.write(g + "\n")
    geneNum += 1
srGenes.close()
