import subprocess
import re

command_line = ['vsearch', '-query',
                "/home/ruben/PycharmProjects/Stage_Wagenigen/Fasta/" + seq + '.fasta',
                '-out', "/home/ruben/PycharmProjects/Stage_Wagenigen/results/" + seq + '_result',
                "-evalue", "0.00001", "-num_alignments", "1", "-task", "megablast", '-num_threads', '4',
                '-db', '/home/ruben/PycharmProjects/Stage_Wagenigen/refdatabase']
subprocess.call(command_line)