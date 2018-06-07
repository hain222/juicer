# Project Juicer

The goal of this program is to evaluate large amounts of genomic sequences, 
(reads primarily) in fasta or fastq format, and group them based off of the 
number of times a given motif that appears at the beginning of a sequence, 
repeats. One way in wish we wish to apply this program is in extracting 
large amounts of telomeric and their adjacent non-telomeric sequences from 
libraries of genomic reads by taking advantage of a repetitive motif 
sequence which identifies telomerws. Additional details as well as 
examples, can be found in the header of the primary source code file.

THIS PROGRAM IS IN EARLY DEVLOPMENT
EVERYTHING IS SUBJECT TO, AND SHOULD BE EXPECTED TO CHANGE

## Getting Started

Download the repository directly from the github webpage, or use git
git clone https://github.com/hain222/juicer

### Prerequisites

This program requires python3 as well as biopython

```
To install python3 on linux
sudo apt-get install python3

To install biopython3 on linux via pip3 (easiest)
sudo apt-get install pip3
sudo pip3 install biopython
```

### Usage

Once the prerequistes are installed, the program will be functional
Run ./juicer.py -h 
To see what the program requires and it's options

## Authors

Harrison Inocencio

## Acknowledgments

Currently Unavailable

