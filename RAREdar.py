"""
Created on Wed May 9 12:08:07 2022

@author: Frank Zhuang
"""

import utility as u

def main():

    MotifList = ['ACTTGA','ACTTGG','ACTTGT','ACTGGA','ACTGGG','ACTGGT','ACTAGA','ACTAGG','ACTAGT']
    ReverseList = ['TGAACT','TGACCT','TGATCT','TGAACC','TGACCC','TGACCA','TGATCT','TGATCC','TGATCA']
    GeneDictionary = u.readFasta('Gene_Sequence_List.fa')
    #print(GeneDictionary)

    hitDictionary,positionDictionary,sequenceDictionary = RAREdar(GeneDictionary, MotifList, ReverseList)#Run RAREdar

    #u.write_output(hitDictionary, 'mockPWM')
    truepositionDictionary = con_coord(positionDictionary)
    #u.write_output(truepositionDictionary, 'mockPWM', data_type = 'position')
    #u.write_output(sequenceDictionary, 'mockPWM', data_type = 'sequence')
    auto_merger(MotifList, ReverseList, hitDictionary, truepositionDictionary, sequenceDictionary)

    return

'''
#RAREdar
#Read along sequences and compare every 6-bp window to every 6-mer motif in the
list of motifs, reports all instances when a 6-bp window matches a predescribed
motif and have a DR 5bps after it
#Input: Dictionary with gene reference number as keys and gene sequence as values
#Input: A list of all 6-mer motifs to be compared against
#Output: Dictionary with gene reference number as keys and number of hits, the
position of every hit, and the sequence of every hit as values.
'''

def RAREdar(dictionary, DRList, RDRList):

    hits = 0 #define object and assign initial value for hits, positions and sequence
    hitDictionary = {}
    positionDictionary = {}
    hitPosition = []
    sequenceDictionary = {}
    hitSequence = []

    DRReverse = [r[::-1] for r in DRList]#Create list of reversed version of RAREs
    RDRReverse = [r[::-1] for r in RDRList]

    for name in dictionary.keys():
        title = name.split()
        sequence = dictionary[name] #assign the sequence of a gene as a string object
        for bp in range(len(sequence)-17):
            window = sequence[bp:bp+6] #create 6-bp window that reads along the sequence string
            repeat = sequence[bp+11:bp+17]
            for dr in range(len(DRList)): #Exhaustive search for every possible RARE occurance on the sequence
                if DRList[dr] == window and window == repeat: #compare every window to every 6-mer in CODING RARE list
                    hits = hits + 1 #record valid hit
                    hitPosition = hitPosition + [bp]
                    hitSequence = hitSequence + [sequence[bp:bp+17]]
                elif RDRList[dr] == window and window == repeat: #compare every window to every 6-mer in COMPLEMENTARY RARE list
                    hits = hits + 1 #record valid hit
                    hitPosition = hitPosition + [bp]
                    hitSequence = hitSequence + [sequence[bp:bp+17]]
                elif DRReverse[dr] == window and window == repeat: #compare every window to every 6-mer in REVERSED RARE list
                    hits = hits + 1 #record valid hit
                    hitPosition = hitPosition + [bp]
                    hitSequence = hitSequence + [sequence[bp:bp+17]]
                elif RDRReverse[dr] == window and window == repeat: #compare every window to every 6-mer in REVERSED COMPLEMENTARY RARE list
                    hits = hits + 1 #record valid hit
                    hitPosition = hitPosition + [bp]
                    hitSequence = hitSequence + [sequence[bp:bp+17]]
        hitDictionary[name] = hits #Assign hit info as values to the gene reference number as keys
        #print(hitDictionary)
        positionDictionary[name] = hitPosition
        sequenceDictionary[name] = hitSequence
        hits = 0 #reinitialize all objects
        hitPosition = []
        hitSequence = []

    return hitDictionary, positionDictionary, sequenceDictionary

'''
#CONCOORD
#Takes relative space coordinates of found RAREs and convert them to true coordinates
#Input :: a dictionary of gene name to relative coordinates
#Output :: a dictionary of gene name to true coordinates
'''

def con_coord(dictionary):
    newdict = {}
    for name in dictionary.keys():
        title = name.split() #Add ID info to every RARE hit
        start = int(title[3])
        end = int(title[4])
        status = int(title[5])
        truecoord = []
        if status == 1: #Situation for gene on forward strand
            for coord in dictionary[name]:
                truecoord = truecoord + [start+coord] #True coordinate by adding relative coordinate to gene starting coordinate
        elif status == -1: #Situation for gene on reversed starnd
            for coord in dictionary[name]:
                truecoord = truecoord + [end-coord-17] #True coordinate by subtracting relative coordinate to gene ending coordinate
        newdict[name] = truecoord
    return newdict

'''
#AUTO MERGER
#Takes all output from RAREdar and create a single, tab-delimited file of all data types
#Input :: 3 dictionaries, representing hit count, position and sequence
#Output :: tab-delimited plain text file with position and sequence data
'''

def auto_merger(DRList, RDRList, hitDictionary, positionDictionary, sequenceDictionary):

    DRReverse = [r[::-1] for r in DRList]
    RDRReverse = [r[::-1] for r in RDRList]

    fileName = 'RAREdar_Results.txt'
    fileWrite = open(fileName,'w')
    fileWrite.write('Chromosome'+'\t'+'Gene'+'\t'+'Mode'+'\t'+'Coordinate'+'\t'+'Sequence'+'\n')
    Entry = ''
    for name in hitDictionary.keys():
        title = name.split()
        chromosome = title[2]
        geneName = title[1]
        coordList = positionDictionary[name]
        sequenceList = sequenceDictionary[name]
        for i in range(len(coordList)):
            Entry = ''
            if any(sequenceList[i][0:6] == dr for dr in DRList):
                mode = 'Forward Coding'
            elif any(sequenceList[i][0:6] == dr for dr in RDRList):
                mode = 'Forward Complement'
            elif any(sequenceList[i][0:6] == dr for dr in DRReverse):
                mode = 'Reversed Coding'
            elif any(sequenceList[i][0:6] == dr for dr in RDRReverse):
                mode = 'Reversed Complement'
            else:
                mode = 'Exception'
            Entry = Entry +'\t' + chromosome +'\t' + geneName +'\t' + mode +'\t' + str(coordList[i]) +'\t' + sequenceList[i] + '\n'
            fileWrite.write(Entry)
    return

'''
#DR SLIDER
#Use a target pattern and a gap width to find all matching hit of this target within a sequence
#Input :: a target pattern, a genetic sequence, the length of gap between direct repeats
#Output :: number of hits, list of position of all hits, list of sequence of all hits
'''

def dr_slider(target, sequence, gap):
    hits = 0
    hitSequence = []
    hitPosition = []
    targetLength = len(target)*2+gap
    for bp in range(len(sequence)-targetLength):
        window = sequence[bp:bp+len(target)] #create window that reads along the sequence string
        repeat = sequence[bp+targetLength-len(target):bp+targetLength]
        if window == repeat and window == target:
            hits = hits + 1 #record valid hit
            hitPosition = hitPosition + [bp]
            hitSequence = hitSequence + [sequence[bp:bp+17]]
    return hits, hitPosition, hitSequence

main()
