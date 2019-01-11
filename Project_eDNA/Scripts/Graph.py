from collections import defaultdict
import matplotlib.pyplot as plt
import re

file = open("Names.txt", "r")
Taxa = file.readlines()
d = defaultdict(int)
acession = re.compile("[A-Z]{2}_?[0-9]{6}\.1")

for item in Taxa:
    for word in item.split("\n"):
        d[word] += 1

del d['']
plt.ylabel("Frequency")
plt.xlabel("Organism")
plt.bar(range(len(d)), list(d.values()), align="center")
plt.xticks(range(len(d)), list(d.keys()), rotation="vertical")
plt.show()