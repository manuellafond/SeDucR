import re

counter = 0


def countValue():
    global counter
    counter += 1
    return counter


def strtr(s, repl):
    pattern = "|".join(map(re.escape, sorted(repl, key=len, reverse=True)))
    return re.sub(pattern, lambda m: f"{repl[m.group()]}_{countValue()}", s)


file = open("rt_segdup.txt")
input = file.readline()
file.close()

# parse input for segrec
srSpecies = open("rt-sp.txt", "w")
srSpecies.write(re.sub(r"\w+", r"'\g<0>'", input[4 : input.find('"', 5)] + ";"))
srSpecies.close()

srGenes = open("rt-genes.txt", "w")
genesList = input[input.find('"', 5) + 1 :].split(" -G ")[1:-1]
for gtree in genesList:
    maps = gtree[gtree.find('"', 1) + 3 : -1]
    maps = dict(e.split(":") for e in maps.split(" "))

    gtree = gtree[1 : gtree.find('"', 1)]
    gtree = strtr(gtree, maps)

    srGenes.write(f"{gtree};\n")
srGenes.close()
