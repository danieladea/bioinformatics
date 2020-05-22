from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))

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


"""
    TODO: Use this space to implement any additional functions you might need

"""

#take readpairs and put each as a separate kmer in our list
def parseKmers(kmers, readPair, k):
    for x in range(2):
        read = readPair[x]
        for y in range(len(read)-k+1):
            if read[y:y+k] in kmers.keys():
                kmers[read[y:y+k]] += 1
            else:
                kmers[read[y:y+k]] = 1
    return

#create list of possible path starts by using kmers that appear often enough to not be errors
def getPaths(reads, minimum, k):
    kmers={}
    allPaths={}
    for readPair in reads:
        parseKmers(kmers, readPair, k)
    for kmer in kmers:
        if kmers[kmer] > minimum:
            if (kmer[:-1] in allPaths.keys()) and (kmer[1:] not in allPaths[kmer[:-1]]):
                allPaths[kmer[:-1]].append(kmer[1:])
            else:
                allPaths[kmer[:-1]] = [kmer[1:]]
    print("Got Paths")
    return allPaths

# check each path to create contigs
def assembler(allPaths):
    #use in and out degrees to get nodes that don't branch
    outdegrees={}
    indegrees={}
    for key,vals in allPaths.items():
        if key not in outdegrees:
            outdegrees[key] = len(vals)
        else:
            outdegrees[key]+=len(vals)
        for val in vals:
            if val not in indegrees:
                indegrees[val] = 1
            else:
                indegrees[val]+=1            
    nodes = list(set(list(indegrees.keys()) +list(outdegrees.keys())))
    nonbranchNodes = []
    branchNodes = []
    for node in allPaths:  
        if (node in indegrees.keys()) and (node in outdegrees.keys()):
            if indegrees[node] ==1 and outdegrees[node]==1:
                nonbranchNodes.append(node)
            else:
                branchNodes.append(node)
        else:
            branchNodes.append(node)

    #use helper function to get contigs 
    contList=[]
    tempPaths = allPaths.copy()
    for path in branchNodes:
        getContig(nonbranchNodes, tempPaths, contList, path)

    #place contigs in a list for output helper function
    contigs = []
    for contig in contList:
        final = ""
        for i in range(len(contig)-1):
            final += contig[i][0]
        final += contig[-1]
        contigs.append(final)

    print("Got contigs")   
    return contigs

#get a contig by using algorithm that finds sequences of nonbranching nodes   
def getContig(nonbranchNodes, allPaths, contList, start):
    
    #in case it's checking one that's not in our adjlist
    #check if we can still find a contig
    try:
        numContigs = len(allPaths[start])
    except KeyError:
        numContigs=0

    if numContigs > 0:
        #find each contig
        for i in range(numContigs):
            contig = []
            contig.append(start)
            nextNode = allPaths[start][0]
            contig.append(nextNode)
            allPaths[start].remove(nextNode)
            if(not allPaths[start]):
                del allPaths[start]

            while nextNode in nonbranchNodes:
                nonbranchNodes.remove(nextNode)
                if nextNode in allPaths.keys():
                    temp = nextNode
                    nextNode = allPaths[nextNode][0]
                    allPaths[temp].remove(nextNode)
                    if(not allPaths[temp]):
                        del allPaths[temp]
                    contig.append(nextNode)
                else:
                    break
            
            contList.append(contig)
            #print(+1 contig)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file
    #24-mers with minimum appearances of kmer = 3 to factor for errors
    k = 24
    minimum=3
    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """
    allPaths = getPaths(input_reads, minimum, k)
    contigs=assembler(allPaths)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
