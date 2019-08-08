import random
import os
#making the sequence alignments folder
if not os.path.exists("sequence_alignments/"):
    os.makedirs("sequence_alignments")
#getting the list of seed sequences
SEQ_LIST = []
SEQ_DICT = {}
for fileName in ["../hg19CCANew_intron_remove2.txt", "../mitochondria_sequences.txt"]:
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
fileNameList = []
for s in ["1", "2"]:
    fileNameList.append("RAW_DATA/HEK293T_DM_" + s + ".fasta")

toAlignLater = {}
class Sequence:
    def __init__(self, seq, ID):
        self.seq = seq
        self.ID = ID
    def getSeq(self):
        return self.seq
    def getID(self):
        return self.ID

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
            if mismatch > 3:
                return False
    return True

for fileName in fileNameList:
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
                    else:
                        if checkSeq(seq, SEQ[:-3] + "CTTTGAGCCTAATGCCTGAA" + SEQ[-3:]):
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
    output = open("sequence_alignments/" + SEQ + ".txt", "w")
    for elem in SEQ_DICT[SEQ]:
        output.write(elem.getID() + elem.getSeq() + "\n")
    output.close()
