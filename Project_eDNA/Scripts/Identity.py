import re
import os
import shutil


list = []


def main():
    hitter()
    listmk()


def hitter():
    src_dict = "/home/ruben/PycharmProjects/Stage_Wagenigen/results/"  # Specify base directory
    result_dict = "/home/ruben/PycharmProjects/Stage_Wagenigen/90%/"

    pattern = re.compile("\(9\d%\)")  #CPatter to search for 90% and higher
    Name_organism = re.compile(">\s[A-Z]{1,}_?[0-9]{5,}\.[0-9]\s[A-Z][a-z]{2,}\s[a-z]{2,}")  #CPatter to search for Name
    Homo_sapiens = re.compile(">\s[A-Z]{1,}_?[0-9]{5,}\.[0-9]\sHomo\ssapiens")  #CPattern to search for acession and Homo sapiens

    for filename in os.listdir(src_dict):
        path = os.path.join(src_dict, filename)
        x = open(path, "r")
        y = open(path, "r")
        z = open(path, "r")

        regex = re.findall(pattern, x.read())
        regex1 = re.findall(Name_organism, y.read())
        regex2 = re.findall(Homo_sapiens, z.read())
        x.close()
        y.close()
        z.close()

        if not regex:
            print("No hit!")
        else:
            print("Hit!")
            if not regex2:
                shutil.move(path, result_dict)
                print("moved", path)
                list.append(regex1)


def listmk():


    file = open("Names.txt", "w")
    with open('Names.txt', 'a') as f:
        for item in list:
                print(item)
                f.write(str(item).replace("]","").replace("[","").replace("'","").replace(">", ""))
                f.write("\n")


main()

