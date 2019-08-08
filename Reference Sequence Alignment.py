## Aligns the list of reference sequences to the human tRNA data set.
## Start Date: 8/5/2019
## Late Edited: 8/8/2019
import random
import os
#making the sequence alignments folder
if not os.path.exists("sequence_alignments/"):
    os.makedirs("sequence_alignments")

class Sequence:
    def __init__(self, seq, ID):
        self.seq = seq
        self.ID = ID
    def getSeq(self):
        return self.seq
    def getID(self):
        return self.ID


#Setting up variables
fileNameList = ["hg19CCA+mito.fasta"]
keepZeroes = True
RAW_fileNameList = []
for f in ["control", "treated"]:
    for s in ["1", "2"]:
        RAW_fileNameList.append("RAW_DATA/tRNA_human_" + f + "_" + s + ".fasta")
#getting the list of seed sequences
SEQ_LIST = []
SEQ_DICT = {}
for fileName in fileNameList:
    handle = open(fileName, "r")
    for line in handle:
        if line[0] != ">":
            if line[:-1] != "":
                SEQ_LIST.append(line[:-1])
                SEQ_DICT[line[:-1]] = []
    handle.close()

UNIQUE = []
for elem in SEQ_LIST:
    if elem in UNIQUE:
        pass
    else:
        UNIQUE.append(elem)
SEQ_LIST = UNIQUE
print(len(SEQ_LIST))

toAlignLater = {}

def checkSeq(seq1, seq2):
    if seq1[-1] == "C":
        seq2 = seq2[:-1]
    if len(seq1) > len(seq2):
        seq1 = seq1[len(seq1) - len(seq2):]
    else:
        seq2 = seq2[len(seq2) - len(seq1):]
    mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
            if mismatch > 2:
                return False
    return True

for fileName in RAW_fileNameList:
    print(fileName)
    i = 0
    handle = open(fileName, "r")
    IDline = ""
    for line in handle:
        if line[0] == ">":
            i += 1
            if i % 10000 == 0:
                print(i/1000)
            IDline = line
        else:
            seq = line[:-1]
            if seq != "":
                matches = []
                for SEQ in SEQ_LIST:
                    if checkSeq(seq, SEQ):
                        matches.append(SEQ)
                if len(matches) == 1:
                    SEQ_DICT[matches[0]].append(Sequence(seq, IDline))
                elif len(matches) > 1:
                    stringMatches = ",".join(matches)
                    if stringMatches in toAlignLater:
                        toAlignLater[stringMatches].append(Sequence(seq, IDline))
                    else:
                        toAlignLater[stringMatches] = [Sequence(seq, IDline)]
    handle.close()
print(len(toAlignLater))
for stringmatches in toAlignLater:
    matches = stringmatches.split(",")
    print(matches)
    seq_list = toAlignLater[stringmatches]
    numSeqeunces = len(seq_list)
    counts = []
    given = []
    for elem in matches:
        counts.append(int(len(SEQ_DICT[elem])))
        given.append(0)
    num_seq = 0
    if sum(counts) > 0:
        for i in range(len(matches)):
            num = int(float(counts[i])/float(sum(counts)) * float(numSeqeunces))
            for it in range(num):
                SEQ_DICT[matches[i]].append(seq_list[num_seq])
                given[i] += 1
                num_seq += 1
    if num_seq != len(seq_list):
        for i in range(len(seq_list) - num_seq):
            biggestError = 0
            toAdd = ""
            if sum(counts) > 0 and sum(given) > 0:
                for j in range(len(matches)):
                    error = float(counts[j])/float(sum(counts)) - float(given[j])/float(sum(given))
                    if error > biggestError:
                        toAdd = matches[j]
                        biggestError = error
            if toAdd == "":
                toAdd = matches[random.randint(0, len(matches) - 1)]
            SEQ_DICT[toAdd].append(seq_list[num_seq])
            num_seq += 1

for SEQ in SEQ_LIST:
    keep = True
    if not keepZeroes:
        if len(SEQ_DICT[SEQ]) == 0:
            keep = False
    if keep:
        output = open("sequence_alignments/" + SEQ + ".txt", "w")
        for elem in SEQ_DICT[SEQ]:
            output.write(elem.getID() + elem.getSeq() + "\n")
        output.close()
        
