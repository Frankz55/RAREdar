# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:06:29 2020

@author: Frank
"""


'''
#BP READER
#Read the extracted kmer and reference pattern side by side to assess match
#Input :: Extracted kmer, reference pattern, arguement of existance of gaps in
#the reference pattern
#Output :: Arguement of match between the extracted kmer and reference pattern
'''
def bpReader(kmer, pattern, gaps): #function to read two sequences side by side
    hit = 0 #reset hit count

    if len(kmer) == len(pattern): #check length of the two sequences
        for i in range(len(kmer)): #for every base in sequence 1
            if kmer[i] == pattern[i]: #if base match for sequence 1 and 2
                hit = hit + 1 #add 1 to hit count
            if gaps == True: #if gaps exist
                if pattern[i] == '-': #if sequence 2 is a gap
                    hit = hit + 1 #add 1 to hit count

    if hit == len(kmer): #if every base matches between the two sequences
        match = True #set match to true
    else:
        match = False #otherwise set match to false

    return match #return boolean arguement



'''
#SEQUENCE GENERATOR
#Generate the pattern with DRs and spaces as a reference for finding the
#promotor in the target region of genes
#Input :: A string used as a pattern, width of the gap, number of DRs
#Output :: Reference pattern as a string
'''
def sequenceGenerator(pattern = 'aggtca', space = 5, repeats = 2): #generate model sequence
    sequence = '' #set sequence to a blank string

    for i in range(repeats): #for every repeat designated
        sequence = sequence + pattern #add pattern to sequence
        if i != repeats-1: #if not the last repeat cycle
            for j in range(space): #for gap length designated
                sequence = sequence + '-' #add gap symbol to sequence

    return sequence



'''
#READ FASTA
#Reads a txt/fasta file of gene names and target regions and compile the data
#Input :: File name
#Output :: Dictionary of gene name to target region
'''
def readFasta(file):
    dictionary = {} #create empty dictionary
    geneNumber = 0 #count the number of genes

    fileHandle = open(file) #open file for reading
    fileLining = fileHandle.readlines() #assemble every line into a list of lines

    sequence = ''

    for line in fileLining: #for every line in the fileLining list
        if line[0] == '>': #if line indicates a gene name
            if sequence != '':
                dictionary[name] = sequence #add name to key of dictionary
                geneNumber = geneNumber + 1 #add 1 to gene number
            name = line.strip() #strip line for name and add to key of dictionary
            sequence = '' #create empty value under key
        else:
            sequence = sequence + line.strip() #strip and add gene sequence line to the new value

    if sequence != '':
        dictionary[name] = sequence
        geneNumber = geneNumber + 1

    return dictionary



'''
#WRITE OUTPUT
#Convert a dictionary output to a space-delimited txt file
#Input :: dictionary of gene name to number of hits, reference pattern, gap
#width, arguement to include genes with 0 hits
#Output :: a txt file pf gene name to number of hits
'''
def write_output(dictionary,pattern,data_type = 'hits',space = 5,zero_check = True,):#write hardFinder results to a file
    fileName = 'pattern_'+str(pattern)+'_spaces_'+str(space)+'_'+data_type+'.txt' #creation of file name with variable pattern and spaces
    fileWrite = open(fileName,'w') #create new file for writing
    fileWrite.write('Genome'+'\t'+data_type+'\n')
    if zero_check == True: #if neglecting promotors with 0 hits
        for items in dictionary.keys(): #for every key in dictionary
            if dictionary[items] != 0 and dictionary[items] != []: #if value of the key is not 0 (0 hits)
                fileWrite.write(items+'\t'+str(dictionary[items])+'\n') #write gene name and hit count to file
    else: #if NOT neglecting promotors with 0 hits
        for items in dictionary.keys(): #for every key in dictionary
            fileWrite.write(items+'\t'+str(dictionary[items])+'\n') #write gene name and hit count to file
    fileWrite.close() #close the new file
    return

'''
#REVERSE COMPLEMENT
#Convert DNA sequences in a dictionary to its reverse complementary strand, 5'->3'
#Input :: dictionary of gene name to sequence
#Output :: dictionary of gene name to sequence, with reverse complementary strand added
'''

def reverse_complement(dictionary):
    newdict = {}
    for name in dictionary.keys():
        sequence = dictionary[name]
        reversed = ''
        basepair = len(sequence)-1
        while basepair != 0:
            base = sequence[basepair]
            if base == 'A':
                reversed = reversed + 'T'
            elif base == 'T':
                reversed = reversed + 'A'
            elif base == 'G':
                reversed = reversed + 'C'
            else:
                reversed = reversed + 'G'
            basepair = basepair - 1
        newname = name + ' Rev'
        newdict[name] = sequence
        newdict[newname] = reversed
    return newdict


def tpc_write_output(dictionary,pattern,space = 5):
    fileName = 'pattern_'+str(pattern)+'_spaces_'+str(space)+'_3pc.txt'
    fileWrite = open(fileName,'w')
    fileWrite.write('Prefix'+'\t'+'Middle'+'\t'+'Suffix'+'\n')
    for index in dictionary.keys():
        if dictionary[index] != []:
            for items in range(len(dictionary[index])):
                fileWrite.write(dictionary[index][items][0:6]+'\t'+dictionary[index][items][6:6+space]+'\t'+dictionary[index][items][6+space:12+space]+'\n')
    fileWrite.close()
    return

def reversal(sequence):
    return sequence[::-1]

def printList(item):
    for i in range(len(item)):
        print(item[i])
    return
