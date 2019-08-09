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
RAW_fileNameList.remove("RAW_DATA/tRNA_human_control_1.fasta")

#getting the list of seed sequences
SEQ_LIST = []
SEQ_DICT = {}
for fileName in fileNameList:
    handle = open(fileName, "r")
    for line in handle:
        if line[0] != ">":
            seq = line[:-1]
            if "C" in seq:
                while seq[-1] != "A" and seq[-1] != "C":
                    seq = seq[:-1]
                SEQ_LIST.append(seq)
                SEQ_DICT[seq] = []
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
##        print("here")
##        print(seq1)
##        print(seq2)
##        print(seq2[-1])
        seq2 = seq2[:-1]
##        print(seq2)
    if len(seq1) > len(seq2):
        seq1 = seq1[len(seq1) - len(seq2):]
    else:
        seq2 = seq2[len(seq2) - len(seq1):]
    mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
            if mismatch > 3:
                return (False, 0)
    return (True, mismatch)

checked = 0
accepted = 0
for fileName in RAW_fileNameList:
    print(fileName)
    i = 0
    handle = open(fileName, "r")
    IDline = ""
    for line in handle:
        if line[0] == ">":
            if checked % 10000 == 0:
                print(checked, accepted)
##            i += 1
##            if i % 10000 == 0:
##                print(i/1000)
            IDline = line
        else:
            seq = line[:-1]
            if seq != "":
                matches = []
                checked += 1
                if "C" in seq:
                    while seq[-1] != "A" and seq[-1] != "C":
                        seq = seq[:-1]
                    lowest_misses = 4
                    for SEQ in SEQ_LIST:
                        check_output = checkSeq(seq, SEQ)
                        if check_output[0]:
                            if check_output[1] < lowest_misses:
                                matches = []
                                matches.append(SEQ)
                                lowest_misses = check_output[1]
                            elif check_output[1] == lowest_misses:
                                matches.append(SEQ)
                    if len(matches) == 1:
                        accepted += 1
                        SEQ_DICT[matches[0]].append(Sequence(seq, IDline))
                    elif len(matches) > 1:
                        accepted += 1
                        stringMatches = ",".join(matches)
                        if stringMatches in toAlignLater:
                            toAlignLater[stringMatches].append(Sequence(seq, IDline))
                        else:
##                            print(len(toAlignLater))
                            toAlignLater[stringMatches] = [Sequence(seq, IDline)]
    handle.close()
print(len(toAlignLater))
for stringmatches in toAlignLater:
    matches = stringmatches.split(",")
##    print(matches)
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
output = ""
for SEQ in SEQ_LIST:
    keep = True
    if not keepZeroes:
        if len(SEQ_DICT[SEQ]) == 0:
            keep = False
    if keep:
        if os.path.exists("sequence_alignments/" + SEQ + ".txt"):
            output = open("sequence_alignments/" + SEQ + ".txt", "a")
            print("Appending...")
        else:
            output = open("sequence_alignments/" + SEQ + ".txt", "w")
        for elem in SEQ_DICT[SEQ]:
            output.write(elem.getID() + elem.getSeq() + "\n")
        output.close()
        
