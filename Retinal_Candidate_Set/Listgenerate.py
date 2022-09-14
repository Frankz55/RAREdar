

def main():

    fileHandle = open('D_rerio_RetinalGene.csv')
    fileLining = fileHandle.readlines()

    geneList = []

    for line in fileLining:
        line.strip()
        geneLine = line.split(',')
        for gene in geneLine:
            if gene not in geneList:
                geneList.append(gene)

    fileName = 'D_rerio_Retinal_Set.txt'
    fileWrite = open(fileName,'w')

    for gene in geneList:
        fileWrite.write(gene+'\n')

    return

main()
