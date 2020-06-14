import sys
from Bio import SeqIO

sequence_to_retrieve = set()

output = open(sys.argv[1] + ".fasta", "w")

for line in open(sys.argv[1]):
	sequence_to_retrieve.add(line.strip("\n"))

for record in SeqIO.parse(open(sys.argv[2]), "fasta"):
	if record.id in sequence_to_retrieve:
		output.write(">" + record.description + "\n" + str(record.seq) + "\n")

