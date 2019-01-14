# Project envoirmental DNA

This is my biodiversity analyser tool, to detect fish in the North Sea. With the use of envoirmental DNA samples

Programs needes:
- VSEARCH
  git clone https://github.com/torognes/vsearch.git
  cd vsearch
  ./autogen.sh
  ./configure
  make
  
- Porechop
  git clone https://github.com/rrwick/Porechop.git
  cd Porechop
  python3 setup.py install

- Filtlong
  git clone https://github.com/rrwick/Filtlong.git
  cd Filtlong
  make -j

- Local BLAST DB

- Mathplotlib
  python3 -mpip install matplotlib
- Biopython
  pip3 install biopython
- Pyqt5
  pip3 install pyqt5


>This Readme is still under construction
