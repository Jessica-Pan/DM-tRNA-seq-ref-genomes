#Creates the outputs specified in the word document tRNA-seq output 080719
#starts at step 4
import os
if not os.path.exists("tRNA-seq-outputs/"):
    os.makedirs("tRNA-seq-outputs")
#Getting the list of seed sequences and creating a dictionary for their isodecoders

SEQ_LIST = []
ISO_DICT = {}

print("getting the seed sequence and identifiers")
name_list = ["SRR1836123", "SRR1836124", "SRR1836125", "SRR1836126"]
short_names = {"SRR1836123": "hc1", "SRR1836126": "ht2", "SRR1836124": "hc2", "SRR1836125": "ht1"}
fancy_names = {"SRR1836123": "human_control_1", "SRR1836126": "human_treated_2", "SRR1836124": "human_control_2", "SRR1836125": "human_treated_1"}
indexDict = {"SRR1836123": 0, "SRR1836126": 3, "SRR1836124": 1, "SRR1836125": 2}
nucIndexes = {"A": 0, "C": 1, "G": 2, "T": 3, "N":4}

handle = open("hg19CCA+mito.fasta", "r")
ISO = ""
for line in handle:
    if line[0] == ">":
        if line[1:3] == "mt" or len(line.split()) == 1:
            ISO = line[1:-1]
        else:
            ISO = line[line.find("-") + 1: line.find("-") + 7] + "_c" + \
                  line[line.find("chr") + 3:line.find(".")] + "t" + \
                  line[line.find("trna") + 4:line.find("-")]
    else:
        if line[0] != "\t":
            seq = line[:-1]
            while seq[-1] != "C" and seq[-1] != "A":
                seq = seq[:-1]
            SEQ_LIST.append(seq)
            ISO_DICT[seq] = ISO
handle.close()
total_raw_counts = [0] * len(fancy_names)
print("getting total counts in raw data")
#getting the total counts per file in the original raw data
for fileName in os.listdir("RAW_DATA"):
    if fileName.endswith(".fasta"):
        handle = open("RAW_DATA/" + fileName, "r")
        count = 0
        done = False
        for line in handle:
            if line[0] == ">":
                if not done:
                    done = True
                    index = indexDict[line[1:line.find(".")]]
                total_raw_counts[index] += 1
        handle.close()

print("going through sequence_alignments")
#going through the data in the sequence_alignments folder in order to create outputs #1, 2, and data set quality
COUNT = 0
IsodecoderCounts = {} # dicts of decoder: [file_1_count, file2_count, ...]
acceptorCounts = {}
alignedCounts = [0] * len(fancy_names)

for fileName in os.listdir("sequence_alignments"):
    if fileName.endswith(".txt"):
        COUNT += 1
        if COUNT % 100 == 0:
            print(COUNT)
        SEQ = fileName[:-4]
        handle = open("sequence_alignments/" + fileName)
        fileCounts = [0] * len(fancy_names)
        for line in handle:
            if line[0] == ">":
                ID = line[1:line.find(".")]
                fileCounts[indexDict[ID]] += 1
        for i in range(len(fileCounts)):
            alignedCounts[i] += fileCounts[i] 
        iso = ISO_DICT[SEQ]
        if iso in IsodecoderCounts:
            for i in range(len(fileCounts)):
                IsodecoderCounts[iso][i] += fileCounts[i]
        else:
            IsodecoderCounts[iso] = fileCounts
            
        anti = iso[:6]
        if anti in acceptorCounts:
            for i in range(len(fileCounts)):
                acceptorCounts[anti][i] += fileCounts[i]
        else:
            acceptorCounts[anti] = fileCounts

        handle.close()

#making the three outputs
print("making outputs")
#Data quality output
output = open("tRNA-seq-outputs/data_quality_output.txt", "w")
output.write("Sample Name\tID\tNumber of total reads\tNumber of reads that aligned to a reference sequence\t% reads aligned\n")
for name in name_list:
    index = indexDict[name]
    percentAligned = float(alignedCounts[index])/float(total_raw_counts[index])
    output.write(name + "\t" + fancy_names[name] + "\t" + str(total_raw_counts[index]) + \
                 "\t" + str(alignedCounts[index]) + "\t" + str(percentAligned) + "\n")
output.close()

#Isodecoder output
#At the same time, making the list of isodecoders that fulfill the requirements to create a mutation fraction per position output
mutationIsodecoders = []
iso_total = [0] * len(name_list)
for j in range(len(iso_total)):
    for i in IsodecoderCounts:
        iso_total[j] += IsodecoderCounts[i][j]
output = open("tRNA-seq-outputs/isodecoder_fraction_output.txt", "w")
output.write("Isodecoder\t" + "\t".join([fancy_names[i] for i in name_list]) + "\n")
decoderFractions = [0] * len(name_list)
biggestDecoderFractions = [0] * len(name_list)
for decoder in IsodecoderCounts:
    for i in range(len(name_list)):
        decoderFractions[i] = float(IsodecoderCounts[decoder][i])/float(iso_total[i])
        if decoderFractions[i] > biggestDecoderFractions[i]:
            biggestDecoderFractions[i] = decoderFractions[i]
    output.write(decoder + "\t" + "\t".join([str(x) for x in decoderFractions]) + "\n")
for decoder in IsodecoderCounts:
    for i in range(len(name_list)):
        if IsodecoderCounts[decoder][i] > 50:
            if float(IsodecoderCounts[decoder][i])/float(iso_total[i]) >0.001 * float(biggestDecoderFractions[i]):
                if decoder not in mutationIsodecoders:
                    mutationIsodecoders.append(decoder)
                break
output.close()

#Acceptor output
acc_total = [0] * len(name_list)
for j in range(len(acc_total)):
    for i in acceptorCounts:
        acc_total[j] += acceptorCounts[i][j]
output = open("tRNA-seq-outputs/anticodon_fraction_output.txt", "w")
output.write("Anticodon\t" + "\t".join([fancy_names[i] for i in name_list]) + "\n")
acceptorFractions = [0] * len(name_list)
for anticodon in acceptorCounts:
    for i in range(len(name_list)):
        acceptorFractions[i] = str(float(acceptorCounts[anticodon][i])/float(acc_total[i]))
    output.write(anticodon + "\t" + "\t".join(acceptorFractions) + "\n")
output.close()

#Mutations per position output
#At the same time, making the charging outputs
print("making per sequence outputs")
print(len(mutationIsodecoders))
counter = 0
charging_output = open("tRNA-seq-outputs/charging_output.txt", "w")
charing_output.write("Isodecoder\t")
for j in range(len(name_list)):
    charging_output.write("3'CCA count\t3'CC count\tTotal Count\tCharging level\t")
charging_output.write("\n")
for decoder in mutationIsodecoders:
    counter += 1
    if counter % 50 == 0:
        print(counter)
    SEQ = []
    for seq in ISO_DICT:
        if ISO_DICT[seq] == decoder:
            SEQ.append(seq)
    if len(SEQ) > 1:
        print(SEQ)
        print("This isodecoder matches more than one sequence. I'm not qualified to handle that")
    SEQ = SEQ[0]
    output = open("tRNA-seq-outputs/" + decoder + ".txt", "w")
    output.write("Position\tOriginal Seed Sequence\t")
    for name in name_list:
        short = short_names[name]
        output.write(short + "total count\tA\tC\tG\tT\tMutation Fraction\tStop Fraction\tMutation fraction_WC\tMI index\t\t")
    output.write("\n")
    handle = open("sequence_alignments/" + SEQ + ".txt", "r")
    GRAND = []
    for i in range(len(SEQ)):
        GRAND.append([])
        for j in range(len(name_list)):
            GRAND[i].append([0, 0, 0, 0, 0])
    index = 0
    charging_list = []
    for j in range(len(name_list)):
        charging_list.append([0, 0])
    for line in handle:
        if line[0] == ">":
            index = indexDict[line[1:line.find(".")]]
        else:
            seq = line[:-1].upper()
            while seq[-1] != "C" and seq[-1] != "A":
                seq = seq[:-1]
            if seq[-1] == "C" or seq[-1] == "A":
                charging_list[index][nucIndexes[seq[-1]]] += 1
            else:
                print(seq)
                print("The ending of this sequence was unexpected and will not be accounted for the charging data")
            if len(seq) > len(SEQ):
                seq = seq[len(seq) - len(SEQ):]
            for pos in range(len(seq)):
                GRAND[-(pos + 1)][index][nucIndexes[seq[-(pos + 1)]]] += 1
    charging_output.write(decoder + "\t")
    for elem in charing_list:
        charging_output.write(str(elem[0]) + "\t" + str(elem[1]) + "\t" + str(sum(elem)) + \
                              "\t" + str(float(elem[0])/float(sum(elem))) + "\t")
    charging_output.write("\n")
    for i in range(len(SEQ)):
        output.write(str(i) + "\t" + SEQ[i] + "\t")
        mutLists = GRAND[i]
        for mutList_index in range(len(mutLists)):
            mutList = mutLists[mutList_index]
            output.write(str(sum(mutList)) + "\t" + "\t".join(str(x) for x in mutList[:-1]) + "\t")
            mutFrac = "0"
            if sum(mutList) > 0:
                mutFrac = (str(float(sum(mutList) - mutList[nucIndexes[SEQ[i]]])/float(sum(mutList))))
            output.write(mutFrac + "\t")
            if i != len(SEQ) - 1:
                nextTotal = sum(GRAND[i + 1][mutList_index])
                if nextTotal > 0:
                    stop = 1 - float(sum(mutList))/float(nextTotal)
                    WC_mut = float(sum(mutList) - mutList[nucIndexes[SEQ[i]]])/float(nextTotal)
                    output.write(str(stop) + "\t" + str(WC_mut) + "\t" + str(WC_mut + stop) + "\t\t")
                else:
                    output.write("\t\t\t\t")
            else:
                output.write("\t\t\t\t")
        output.write("\n")
    output.close()

print("done :)")
