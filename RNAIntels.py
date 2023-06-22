#!/bin/env python3
import argparse
import numpy as np
import pandas as pd

from lib.main import RNAIntels

def main():
    parser = argparse.ArgumentParser(
                    prog='RNAIntels',
                    description='Classify RNA sequences into coding and noncoding',
                    epilog='Please cite us')
    parser.add_argument('sequence_file', type=str)
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="Enables the verbose mode")
    parser.add_argument('-m', '--model', type=str, default=None, help="Select the applied .h5 file")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Defines how many CPU threads should be used")
    parser.add_argument('--min', type=int, default=25, help="Minimal sequence length")
    parser.add_argument('--max', type=int, default=60, help="Maximal sequence length")
    args = parser.parse_args()
    RNAIntels(args)

if __name__ == "__main__":
    main()