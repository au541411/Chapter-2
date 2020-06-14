import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

handle = open(sys.argv[1])
protein_seq_file = open(sys.argv[1] + ".faa", "w")

for record in SeqIO.parse(handle, "fasta"):
#	seq_id = record.description.split(" ")[1]
	protein_seq_file.write(">" + record.id + ".0\n" + str(record.seq.translate(table=11)) + "\n")
	protein_seq_file.write(">" + record.id + ".1\n" + str(record.seq[1:].translate(table=11)) + "\n")
	protein_seq_file.write(">" + record.id + ".2\n" + str(record.seq[2:].translate(table=11)) + "\n")
	reversed_seq = record.seq.reverse_complement()
	protein_seq_file.write(">" + record.id + ".0rc\n" + str(reversed_seq.translate(table=11)) + "\n")
	protein_seq_file.write(">" + record.id + ".1rc\n" + str(reversed_seq[1:].translate(table=11)) + "\n")
	protein_seq_file.write(">" + record.id + ".2rc\n" + str(reversed_seq[2:].translate(table=11)) + "\n")

handle.close()
