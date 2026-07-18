import argparse
import os
import statistics

# Parse command-line arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("MSA", help="Path to MSA file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

def defragment_MSA(MSA, output):
    with open(MSA, 'r') as input, open(output, 'w') as output:
        i = 0
        for line in input:
            if (line[0] == ">") and (i==0):
                output.write(line)
                i += 1
            elif line[0] == ">":
                output.write("\n"+line)
            else:
                output.write(line[0:-1])


defragment_MSA(args.MSA, args.output+".def")

input_basename = os.path.basename(args.output)
headers = []
sequences = []
with open(args.output+".def", 'r') as input:
    for line in input:
        if line.startswith('>'):
            headers.append(line)
            sequences.append([])
        else:
            sequences[-1].append(line.strip().upper())

cred = [float(h.split('_')[-1]) for h in headers]
n = len(headers)
cluster_credibility = round(statistics.mean(cred), 2)

sequences = [''.join(seq_lines) for seq_lines in sequences]
sequence_length = len(sequences[0])

consensus = []
# A T G C -
for base in range(sequence_length):
    count = [0, 0, 0, 0, 0]
    for seq in sequences:
        letter = seq[base]
        if letter == "A":
            count[0] += 1
        elif letter == "T":
            count[1] += 1
        elif letter == "G":
            count[2] += 1
        elif letter == "C":
            count[3] += 1
        elif letter == "-":
            count[4] += 1

    m = max(count)
    if m == count[0]:
        consensus.append("A")
    elif m == count[1]:
        consensus.append("T")
    elif m == count[2]:
        consensus.append("G")
    elif m == count[3]:
        consensus.append("C")
    else:
        if count[4] < (count[0]+count[1]+count[2]+count[3]):
            b = [count[0], count[1], count[2], count[3]]
            mb = max(b)
            if mb == count[0]:
                consensus.append("A")
            elif mb == count[1]:
                consensus.append("T")
            elif mb == count[2]:
                consensus.append("G")
            elif mb == count[3]:
                consensus.append("C")
        else:
            consensus.append("-")

with open(args.output, 'w') as output:
    output.write(">" + input_basename + "_" + str(cluster_credibility) + "_" + str(n) + "\n")
    output.write(''.join(consensus))
    output.write("\n")

os.remove(args.output+".def")
