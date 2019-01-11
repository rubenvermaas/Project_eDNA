from Bio import SeqIO

try:
     from StringIO import StringIO
except ImportError:
     from io import StringIO

File = open("Allfastq.fastq","r")
File2 = open("Allfasta.fastq", "w")


handle = StringIO("")
SeqIO.convert(File, "fastq", handle, "fasta")
File2.write(handle.getvalue())

