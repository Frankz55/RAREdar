

def main():

    fileHandle = open('RetinalGene.csv')
    fileLining = fileHandle.readlines()

    geneList = []

    for line in fileLining:
        line.strip()
        geneLine = line.split(',')
        for gene in geneLine:
            if gene not in geneList:
                geneList.append(gene)

    fileName = 'CandidateGenes.txt'
    fileWrite = open(fileName,'w')

    for gene in geneList:
        fileWrite.write(gene+'\n')

    return

main()
