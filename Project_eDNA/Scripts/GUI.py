import sys
import subprocess
from PyQt5.QtWidgets import *
from Bio import SeqIO

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from collections import defaultdict
import matplotlib.pyplot as plt
import re
import shutil
import os
import time
import csv

class window(QMainWindow):

    def __init__(self, parent=None):
        super(window, self).__init__(parent)
        self.setGeometry(550, 250, 850, 550)
        self.setWindowTitle('eDNA Fish Analyzer')

        #Menu bar File and Edit and activates the home function

        goHome = QAction("&Home", self)
        goHome.setShortcut("Ctrl+H")
        goHome.setStatusTip("Go to homepage")
        goHome.triggered.connect(self.home)

        extractAction = QAction('&Quit', self)
        extractAction.setShortcut('Ctrl+Q')
        extractAction.setStatusTip('leave the app')
        extractAction.triggered.connect(self.close_application)

        openEditor = QAction('&Editor', self)
        openEditor.setShortcut('Ctrl+E')
        openEditor.setStatusTip('Open Editor')
        openEditor.triggered.connect(self.editor)

        openFile = QAction('&Open File', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open File')
        openFile.triggered.connect(self.file_open)

        working_dir = QAction('&Working dir', self)
        working_dir.setShortcut('Ctrl+s')
        working_dir.setStatusTip('Select working dir')
        working_dir.triggered.connect(self.working_dir)

        self.statusBar()

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(goHome)
        fileMenu.addAction(openFile)
        fileMenu.addAction(working_dir)
        fileMenu.addAction(extractAction)

        editorMenu = mainMenu.addMenu('&Edit')
        editorMenu.addAction(openEditor)

        self.figure = plt.figure(figsize=(6, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.move(230, 90)
        self.canvas.setParent(self)

        self.home()

    #This function will edit files

    def editor(self):
        self.textEdit = QTextEdit()
        self.setCentralWidget(self.textEdit)

    #This function will save written files

    #def file_save(self):
     #   name, _ = QFileDialog.getSaveFileName(self, 'Save File', options=QFileDialog.DontUseNativeDialog)
     #   file = open(name, 'w')
      #  text = self.textEdit.toPlainText()
      #  file.write(text)
       # file.close()

    #This function will open directory

    def dir_open(self):

        self.fileDir = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.myTextBox.setText(self.fileDir)

        if self.fileDir:
            return self.fileDir

   #This function will open a file

    def file_open(self):
        self.fileName, _ = QFileDialog.getOpenFileName(self, 'Open File', options=QFileDialog.DontUseNativeDialog)
        self.myTextBox.setText(self.fileName)

        if self.fileName:
            return self.fileName

    # This function saves a text file wit the working dir

    def working_dir(self):
        working_dir = open("Working_dir.txt", "w")
        fileDir = QFileDialog.getExistingDirectory(self, 'Select Directory')
        working_dir.write(fileDir)
        working_dir.close()

        working_dir.close()

    #This function controles where all the buttons, checkboxes, labels, etc are located and add them to different functions

    def home(self):

        btn = QPushButton('Load Dir', self)
        btn.clicked.connect(self.dir_open)
        btn.resize(btn.sizeHint())
        btn.move(10, 55)

        btn2 = QPushButton('Load File', self)
        btn2.clicked.connect(self.file_open)
        btn2.resize(btn2.sizeHint())
        btn2.move(10, 25)

        self.btn1 = QPushButton('Run', self)
        self.btn1.resize(btn.sizeHint())
        self.btn1.move(20, 470)
        self.btn1.clicked.connect(self.pipeline)

        self.filtering = QCheckBox("Filtering:", self)
        self.filtering.setChecked(False)
        self.filtering.move(10, 120)

        self.adapters = QCheckBox("Adapters", self)
        self.adapters.setChecked(False)
        self.adapters.move(110, 80)

        self.fastq = QCheckBox("Fastq file", self)
        self.fastq.setChecked(False)
        self.fastq.move(10, 80)

        self.singletons = QCheckBox("Exclude singletons", self)
        self.singletons.setChecked(False)
        self.singletons.resize(150, 25)
        self.singletons.move(10, 255)

        box = self.myTextBox = QTextEdit(self)
        box.resize(732, 25)
        box.move(100,25)

        self.box1 = QLineEdit(self)
        self.box1.resize(50, 25)
        self.box1.move(170, 150)

        self.box2 = QLineEdit(self)
        self.box2.resize(50, 25)
        self.box2.move(170, 180)

        self.box3 = QLineEdit(self)
        self.box3.resize(50, 25)
        self.box3.move(170, 220)

        self.box4 = QLineEdit(self)
        self.box4.resize(100,25)
        self.box4.move(10,320)

        label1 = QLabel(self)
        label1.setText("Sequence divergence")
        label1.resize(150,25)
        label1.move(10, 220)

        label2 = QLabel(self)
        label2.setText("Read length cut-off")
        label2.resize(150, 25)
        label2.move(30, 150)

        label3 = QLabel(self)
        label3.setText("Q-score cut-off")
        label3.resize(150, 25)
        label3.move(30, 180)

        label4 = QLabel(self)
        label4.setText("Database:")
        label4.resize(200, 25)
        label4.move(10, 290)

        self.progress = QProgressBar(self)
        self.progress.setGeometry(227, 500, 275, 20)

    #This function will run the functions in the selected order. Dependent on wich inputs you gave in the GUI. It will also display the progressbar.

    def pipeline(self):

        self.completed = 0

        #open working dir before the analysis is ran
        working_file = open("Working_dir.txt", "r")
        self.working_file1 = working_file.read()

        if self.working_file1 == "":
            self.Pop_up_wkdir()
        else:
            pass

        #pipeline when fastq, filtering, adapter and singletons checkboxes are clicked

        if self.fastq.isChecked() and self.adapters.isChecked() and self.filtering.isChecked() and self.singletons.isChecked():
            while self.completed <100:
                start = time.time()

                #pre-processing
                self.Remove_tmp()
                self.completed += 1.0
                self.progress.setValue(self.completed)
                self.Filtlong()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.Fastq_omzetter()
                self.completed += 4.0
                self.progress.setValue(self.completed)
                self.Porechop()
                self.completed += 25.0
                self.progress.setValue(self.completed)

                #analysis
                self.Clustering()
                self.completed += 25.0
                self.progress.setValue(self.completed)
                self.splitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.Singleton()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.local_Blast()
                self.completed += 20.0
                self.progress.setValue(self.completed)
                self.hitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)

                #results
                self.listmk()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                end = time.time()
                end1 = end/60
                start1 = start/60
                print(end1 - start1)
                self.plot()

        elif self.fastq.isChecked() and self.adapters.isChecked() and self.filtering.isChecked():
            while self.completed <100:
                start = time.time()
                self.Remove_tmp()
                self.completed += 1.0
                self.progress.setValue(self.completed)
                self.Filtlong()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.Fastq_omzetter()
                self.completed += 4.0
                self.progress.setValue(self.completed)
                self.Porechop()
                self.completed += 25.0
                self.progress.setValue(self.completed)
                self.Clustering()
                self.completed += 30.0
                self.progress.setValue(self.completed)
                self.splitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.local_Blast()
                self.completed += 20.0
                self.progress.setValue(self.completed)
                self.hitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.listmk()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                end = time.time()
                end1 = end/60
                start1 = start/60
                print(end1 - start1)
                self.plot()

        elif self.adapters.isChecked() and self.filtering.isChecked() and self.singletons.isChecked():
            while self.completed < 100:
                self.Remove_tmp()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.Porechop()
                self.completed += 30.0
                self.progress.setValue(self.completed)
                self.Clustering()
                self.completed += 30.0
                self.progress.setValue(self.completed)
                self.splitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.local_Blast()
                self.completed += 20.0
                self.progress.setValue(self.completed)
                self.hitter()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.listmk()
                self.completed += 5.0
                self.progress.setValue(self.completed)
                self.plot()

        elif self.checkBox.isChecked():
            self.plot()

        else:
            self.Clustering()
            self.splitter()
            self.local_Blast()
            self.hitter()
            self.listmk()
            self.plot()

        # This function will remove all temporaly files if at the start of a new analisis, so that the resultst will not overlap

    def Remove_tmp(self):  # works

        Barcodes = self.working_file1 + '/Barcodes'
        Fasta = self.working_file1 + '/Fasta/'
        Results = self.working_file1 + '/Results/'
        best_hits = self.working_file1 + '/90%/'
        Singletons = self.working_file1 + '/Singletons/'

        for the_file in os.listdir(Fasta):
            file_path = os.path.join(Fasta, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

        for the_file in os.listdir(Barcodes):
            file_path = os.path.join(Barcodes, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

        for the_file in os.listdir(Results):
            file_path = os.path.join(Results, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

        for the_file in os.listdir(best_hits):
            file_path = os.path.join(best_hits, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

        for the_file in os.listdir(Singletons):
            file_path = os.path.join(Singletons, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)

    # This function will call on filtlong on the command  line to set a threshhold for read lenght and for

    def Filtlong(self):

        result = open('inputReads.fastq', 'w')

        for filename in os.listdir(self.fileDir):
            with open(os.path.join(self.fileDir, filename)) as f:
                content = f.readlines()
                for i in content:
                  result.write(i)

        result.close()

        length = self.box1.text()
        score = self.box2.text()

        os.system('filtlong --min_mean_q ' + score + ' --min_length ' + length + ' inputReads.fastq > output.fastq')


    #This function turns a fastq file to an fasta file

    def Fastq_omzetter(self):

        from io import StringIO

        File = open("output.fastq", "r")
        self.File2 = open("fasta_file.fasta", "w")

        handle = StringIO("")
        SeqIO.convert(File, "fastq", handle, "fasta")
        self.File2.write(handle.getvalue())

    # This function will call on Porchop on the commandline to trimm the adapters from the reads

    def Porechop(self):

        InputPore = "fasta_file.fasta"

        self.command_line = ["porechop", "-i", InputPore,
                             "-b", self.working_file1 + "/Barcodes"]
        subprocess.call(self.command_line)

    # This function will call on VSearch on the commandline to cluster the reads into OTU's

    def Clustering(self):  # works

        ID = str(self.box3.text())

        Input = "output_reads.fasta"  # str(self.fileName)

        self.command_line = ["vsearch", "-cluster_fast", Input,
                            "-id", ID, "-consout", "Test.fasta"]
        subprocess.call(self.command_line)

    # This function splits the consencus clusters in single FASTA files

    def splitter(self):

        i = 0

        with open("Test.fasta", "r") as file:
            data = file.read()
            seqs1 = re.split(">", data)

        for seqs2 in seqs1:
            filename = "result" + str(i)
            result = open(self.working_file1 + "/Fasta/" + filename + ".fasta", "w")
            result.write(">" + seqs2)
            i += 1

    # This function will filter out the clusters with 1 sequence, this option can be selected at the GUI

    def Singleton(self):

        src_dict = self.working_file1 + "/Fasta/"
        result_dict = self.working_file1 +"/Singletons/"

        reseqs = re.compile("(seqs=)[1]\n")  # searches for seq=1

        self.result = seqs = []

        for filename in os.listdir(src_dict):
            path = os.path.join(src_dict, filename)
            x = open(path, "r")
            regex = re.findall(reseqs, x.read())
            if regex:
                shutil.move(path, result_dict)
                # print("moved", path)
            else:
                seqs.append(filename)

            x.close()

    # This function will run a BLAST with the found concensus clusters against the selected local blast database

    def local_Blast(self):

        l = 1

        database = self.box4.text()

        for seq in self.result:
            l += 1
            command_line = ["blastn", "-query", self.working_file1 + "/Fasta/" + seq,
                            "-out", self.working_file1 + "/Results/" + seq,
                            "-evalue", "0.00001", "-num_alignments", "1", "-task", "megablast", "-num_threads", "4",
                            "-db", self.working_file1 + "/Databases/" + database]
            subprocess.call(command_line)


    #This function searches for the results with an identity score of higer then 90. If there are Homo sapiens found that hit for 90% he will filter those out

    def hitter(self):

        self.lijst = List = []

        src_dict = self.working_file1 + "/Results/"  # Specify base directory
        result_dict = self.working_file1 + "/90%/"

        pattern = re.compile("\(9\d%\)")  # CPatter to search for 90% and higher
        Name_organism = re.compile(">\s[A-Z]{1,}_?[0-9]{5,}\.[0-9]\s[A-Z][a-z]{2,}\s[a-z]{2,}")  # Regex to search for Organisme name
        Homo_sapiens = re.compile(">\s[A-Z]{1,}_?[0-9]{5,}\.[0-9]\sHomo\ssapiens")  # Regex to search for acession and Homo sapiens

        for filename in os.listdir(src_dict):
            path = os.path.join(src_dict, filename)
            x = open(path, "r")
            y = open(path, "r")
            z = open(path, "r")

            regex = re.findall(pattern, x.read())
            regex1 = re.findall(Name_organism, y.read())
            regex2 = re.findall(Homo_sapiens, z.read())

            if not regex:
                print("No hit!")
            else:
                print("Hit!")
                if not regex2:
                    shutil.move(path, result_dict)
                    print("moved", path)
                    List.append(regex1)

            x.close()
            y.close()
            z.close()

    #This function makes a list with the organismes

    def listmk(self): #naar cvs file zetten

        file = open("Names.txt", "w")
        with open('Names.txt', 'a') as f:
            for item in self.lijst:
                print(item)
                f.write(str(item).replace("]", "").replace("[", "").replace("'", "").replace(">", ""))
                f.write("\n")

        file2 = open("Names.txt", "r")
        filenew = file2.readlines()

        with open("Results.csv", mode="w") as Results:
            for item2 in filenew:
                Results_writer = csv.writer(Results, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                Results_writer.writerow([item2])

        self.Pop_up()


    #This function will create a plot with the frequency results

    def plot(self):

        file = open("Names.txt", "r")
        Taxa = file.readlines()
        d = defaultdict(int)
        print(d)
        acession = re.compile("[A-Z]{2}_?[0-9]{6}\.1")

        for item in Taxa:
            for word in item.split("\n"):
                d[word] += 1
        print(d.keys())
        del d['']
        ax = self.figure.add_subplot(111)
        ax.set_title("eDNA  classification results")
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Organisme")
        ax.bar(range(len(d)), list(d.values()), align="center")
        ax.set_xticks(range(len(d)))
        ax.set_xticklabels(list(d.keys()), rotation=0)
        self.canvas.draw()

    def style_choice(self, text):
        self.styleChoice.setText(text)
        QApplication.setStyle(QStyleFactory.create(text))


    #This function will give a pop-up message if the analysis is completed

    def Pop_up(self):
        popup = QMessageBox.about(self, "Process", "Analysis is done!")

        return popup


    def Pop_up_wkdir(self):
        popup = QMessageBox.about(self, "Error", "No working directory selected")

        return popup

    #This function will give a pop-up message if you close the application

    def close_application(self):

        choice = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit?", QMessageBox.Yes |
                                     QMessageBox.No, QMessageBox.No)

        if choice == QMessageBox.Yes:
            print('quit application')
            sys.exit()
        else:
            pass


if __name__ == "__main__":  # had to add this otherwise app crashed

    def run():
        app = QApplication(sys.argv)
        Gui = window()
        Gui.show()
        sys.exit(app.exec_())

run()



