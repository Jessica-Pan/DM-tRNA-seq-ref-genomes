#making a better formatted hg19CCANew_intron_remove+mito

handle = open("hg19CCANew_intron_remove+mito.txt", "r")
output = open("hg19CCA+mito.fasta", "w")

for line in handle:
    if len(line.split()) == 2:
        output.write(">" + line.split()[0] + "\n" + line.split()[1] + "\n")
    else:
        output.write(line)
output.close()
handle.close()
