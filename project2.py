import sys
import argparse
import numpy as np
import time
import zipfile
from statistics import mode
from collections import defaultdict
from collections import Counter 


allDict = {}
allInserts = {}
allDels = {}
def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""
#take reads and create kmers out of them 
def splitReads(k, input_reads):
    readsDict={}
    tripleK=3*k
    iterationsRange = 50 - (tripleK)+1
    for readArr in input_reads:
        for x in range(iterationsRange):
            for y in range(2):
                if(readArr[y][x:x+tripleK] in readsDict.keys()):
                    readsDict[readArr[y][x:x+tripleK]] = readsDict[readArr[y][x:x+tripleK]] + 1 
                else:
                    readsDict[readArr[y][x:x+tripleK]] = 1


    for k,v in list(readsDict.items()):
        if v < 3:
           del readsDict[k]
    reads = list(readsDict)
    return reads

#create a lookup table from the reference genome 
def makeDict(k, reference):
    refLen = len(reference)
    table = {}
    for x in range(refLen-k+1):
        kmer = reference[x:x+k]
        if kmer not in table.keys():
            table[kmer]=[x]
            #print(f"new kmer {kmer}")
        else:
            table[kmer].append(x)

    #print(table)
    return table

#takes k-length part of read and compares it to reference to see if it is off by only one substitution
def findSnp(read, position, reference, k):
    errors =0
    #snpList = []
    tempDict ={}
    for currPos in range(len(read)):
        refPos = position+currPos
        currRead = read[currPos]
        currRef = reference[refPos] 
        if(currRead!=currRef):
            errors=errors+1
            tempDict[refPos] = currRead
            #snpList.append(refPos,currRef, currRead)

    if(errors < 2 and errors>0):
        for i in tempDict.keys():
                if i not in allDict.keys():
                    allDict[i]=tempDict[i]

#emulates incrementing two pointers to see if the reference and read differ by 1 insertion
def findInsert(read, position, reference, k):
    tempInserts ={}
    referenceKmer = reference[position:position+k-1]
    refCtr=0
    readCtr=0
    while(True):
        if(refCtr>(k-2)):
            break
        elif referenceKmer[refCtr]== read[readCtr]:
            refCtr=refCtr+1
            readCtr=readCtr+1
        else:
            tempInserts[refCtr+position] = read[readCtr]
            refCtr=refCtr+1 
    if(len(list(tempInserts)) == 1):
        for key in tempInserts.keys():
            if key not in allInserts.keys():
                allInserts[key] = tempInserts[key]

#emulates incrementing two pointers to see if the reference and read differ by 1 deletion
def findDel(read, position, reference, k):
    tempDels ={}
    referenceKmer = reference[position:position+k+1]
    refCtr=0
    readCtr=0
    while(True):
        if(refCtr>(k) or readCtr>k-1):
            break
        elif referenceKmer[refCtr] == read[readCtr]:
            refCtr=refCtr+1
            readCtr=readCtr+1
        else:
            tempDels[refCtr+position] = referenceKmer[refCtr]
            refCtr=refCtr+1 
    if(len(list(tempDels)) == 1):
        for key in tempDels.keys():
            if key not in allDels.keys():
                allDels[key] = tempDels[key]

#takes kmers we got from reads and uses algorithm that checks if 2/3 are perfect matches
#if so, check for snp, insert, or deletion depending on locations of matches
def getSnps(k, allReads, hashTable, reference):
    readLen = len(allReads[0])
    #print(readLen)
    #snps = []
    #inserts = {}

    for read in allReads:
        firstRead = read[0:k]
        secondRead = read[k:2*k]
        thirdRead = read[2*k:3*k]

        #difference between first chunk position and third chunk position in reference should be 10
        #if not, check for insertion or deletion accordingly
        if firstRead in hashTable.keys() and thirdRead in hashTable.keys():
            if ((hashTable[thirdRead][0] - hashTable[firstRead][0]) == 2*k):
                findSnp(secondRead, hashTable[firstRead][0]+k, reference, k)
            if ((hashTable[thirdRead][0] - hashTable[firstRead][0]) == 2*k-1):
                findInsert(secondRead, hashTable[firstRead][0]+k, reference, k)
            if ((hashTable[thirdRead][0] - hashTable[firstRead][0]) == 2*k+1):
                findDel(secondRead, hashTable[firstRead][0]+k, reference, k)
        
        if firstRead in hashTable.keys() and secondRead in hashTable.keys():
            if((hashTable[secondRead][0] - hashTable[firstRead][0]) == k):
                findSnp(thirdRead, hashTable[secondRead][0]+k, reference, k)

        if secondRead in hashTable.keys() and thirdRead in hashTable.keys():
            if((hashTable[thirdRead][0] - hashTable[secondRead][0]) == k): 
                findSnp(firstRead, hashTable[secondRead][0]-k, reference, k)             
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here

    """
    k=10
    allReads = splitReads(k, input_reads)
    #print(allReads)
    hashTable = makeDict(k, reference)
    #print(hashTable)
    getSnps(k, allReads, hashTable, reference)


    snps =[]
    for x in allDict.keys(): 
        snps.append([reference[x], allDict[x], x])
    #snps = [['A', 'G', 3425]]


    insertions =[]
    for x in allInserts.keys(): 
        insertions.append([allInserts[x], x])
    #insertions = [['ACGTA', 12434]]

    deletions =[]
    for x in allDels.keys(): 
        deletions.append([allDels[x], x])
    #deletions = [['CACGG', 12]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
