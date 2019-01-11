import subprocess
import re


seqs = []

def main():
    splitter()
    local_Blast()


def splitter():

    i = 0

    with open("Allfasta.fastq", "r") as file:
        data = file.read()
        seqs1 = re.split(">", data)

    for seqs2 in seqs1:
        filename = "result" + str(i)
        result = open("/home/ruben/PycharmProjects/Stage_Wagenigen/Fasta/" + filename + ".fasta", "w")
        result.write(">" + seqs2)
        seqs.append(filename)
        i += 1

    return(seqs)

def local_Blast():

    l = 0

    for seq in seqs:
        l += 1
        print(l)
        command_line = ['blastn', '-query',
                        "/home/ruben/PycharmProjects/Stage_Wagenigen/Fasta/" + seq + '.fasta',
                        '-out',"/home/ruben/PycharmProjects/Stage_Wagenigen/results/" + seq + '_result',
                        "-evalue", "0.00001", "-num_alignments", "1", "-task", "megablast",  '-num_threads', '4',
                        '-db', '/home/ruben/PycharmProjects/Stage_Wagenigen/refdatabase']
        subprocess.call(command_line)

main()

