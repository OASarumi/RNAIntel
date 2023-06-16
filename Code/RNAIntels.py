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
    parser.add_argument('filepath', type=str)
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-m', '--model', type=str, default=None)
    parser.add_argument('-t', '--threads', type=int, default=1)
    args = parser.parse_args()
    RNAIntels(args)

if __name__ == "__main__":
    main()