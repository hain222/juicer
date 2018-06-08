# Project Juicer

The goal of this program is to evaluate large amounts of genomic sequences, 
(reads primarily) in fasta or fastq format, and group them based off of the 
number of times a given motif, that appears at or near the beginning of a sequence, 
repeats. One way we wish to apply this program is to extract
large amounts of telomeric and their adjacent non-telomeric sequences from 
libraries of genomic reads by taking advantage of a repetitive motif 
sequence which identifies telomeres. Additional details as well as 
examples, can be found in the header of the primary source code file.

THIS PROGRAM IS IN EARLY DEVELOPMENT
EVERYTHING IS SUBJECT TO, AND SHOULD BE EXPECTED TO CHANGE

## Getting Started

Download the repository directly from the github webpage, or use git
```
git clone https://github.com/hain222/juicer
```

### Prerequisites

This program requires python3 as well as biopython.
Strongly recommend installing on a linux virtual machine.

```
# To install python3 on linux
sudo apt-get install python3

# To install biopython on linux via pip3 (easiest)
sudo apt-get install python3-pip
sudo pip3 install biopython
```

### Usage

Once the prerequisites are installed, the program will be functional. 
Run juicer with the help flag to see program requirements and optional 
values
```
./juicer.py -h
```

## Authors

Harrison Inocencio

## Acknowledgments

Currently Unavailable

