import sys
import argparse
import time
import zipfile
import statistics 
from statistics import mode
from collections import defaultdict
from collections import Counter 


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
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
#put read pairs in an array of reads
def splitReads(input_reads):
    reads=[]
    for readArr in input_reads:
        reads.append(readArr[0])
        reads.append(readArr[1])
    return reads

#create a lookup table for k-mers in the reference genome
def makeDict(k, reference):
    refLen = len(reference)
    table = {}
    for x in range(refLen-k+1):
        kmer = reference[x:x+k]
        if kmer not in table.keys():
            table[kmer]=[x]
        else:
            table[kmer].append(x)
    return table

#make a dictionary of potential snps 
#take 50-length reads and split into 5; check if any 1/5th perfectly matches to reference genome
#   if so, find any potential snps in the read
#simplified algorithm, but works due to unrealistically high coverage
def getSnps(k, allReads, hashTable, reference):
    readLen = len(allReads[0])
    #snps = []
    allDict={}
    for read in allReads:
        iterations = readLen/k
        for x in range(iterations):
            currRead = read[x*k:x*k+k]
            #print(currRead)
            if currRead in hashTable.keys():
                #print(currRead)
                #print(reference[hashTable[currRead][0]:hashTable[currRead][0]+10])
                tempDict ={}
                for position in hashTable[currRead]:
                    start = position - (x*k)
                    errors = 0
                    currSnps=[]
                    end = start+50
                    for currNode in range(start,end):
                        try:
                            tempRef = reference[currNode]
                            tempRead = read[currNode-start]
                        except IndexError:
                            break
                        if tempRef!= tempRead:
                            errors = errors+1
                            #tempDict[currNode] = [tempRef,tempRead]
                            tempDict[currNode] = tempRead
                        # else: 
                        #     #currSnps.append([tempRef ,tempRead ,currNode])
                        #     tempDict[currNode] = [tempRef,tempRead]
                        if errors>2:
                            break
                        elif currNode == start+49:
                            #print(currSnps)
                            #snps.extend(currSnps)
                            #allDict.update(tempDict)
                            for i in tempDict.keys():
                                if i in allDict.keys():
                                    allDict[i].append(tempDict[i])
                                else:
                                    allDict[i]=[tempDict[i]]
    return allDict
                             
                 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
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
    allReads = splitReads(input_reads)
    #print(allReads)
    hashTable = makeDict(k, reference)
    #print(hashTable)
    snpsDict = getSnps(k, allReads, hashTable, reference)

    #choose the snp that occurred the most
    snps =[]
    for x in snpsDict.keys():
        
        occurence_count = Counter(snpsDict[x]) 
        bestChoice = occurence_count.most_common(1)[0][0]
        appearances = occurence_count.most_common(1)[0][1] 
        if appearances > 20:
            print(appearances, x,snpsDict[x])
            snps.append([reference[x], bestChoice, x])
    
    #snps = [['A', 'G', 3425]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
