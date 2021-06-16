#!/usr/bin/python3

"""
Author: Fabiana Patalano
Mail: fabiana.patalano97@gmail.com
Last updated: 18/05/2021

"""

from Bio import SeqIO
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from argparse import RawTextHelpFormatter
from Bio.Align import MultipleSeqAlignment


def get_args():
    parser = argparse.ArgumentParser(description=" ".join(__doc__.splitlines()[4:]),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', required=True, dest='msa', help='Input stockholm file')

    parser.add_argument('--output', '-out', '-o', required=True, dest='out', help='output file')
    args = parser.parse_args()

    return args


def convertStocktoa3m(seqA, output):
    f = open(output, "w+")
    for record in SeqIO.parse(seqA, "stockholm"):
         f.write(f">{record.name} \n {record.seq}\n")
    f.close()


if __name__ == "__main__":
    args = get_args()
    print(args)

    seqA = os.path.abspath(args.msa)
    output = os.path.abspath(args.out)

    convertStocktoa3m(seqA, output)
