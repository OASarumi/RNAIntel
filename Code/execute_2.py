#!/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd

from multiprocessing import Semaphore, Manager, Process
from tensorflow.keras.models import load_model
from math import sin, cos, pi, floor, ceil
from Bio import SeqIO


# This FCGR is a variant version of the R package kaos by Dominic Eger and Hannah Franziska Löchel, implemented in Python by Sandra Clemens
class CGR:
    r = 1
    
    #preciion 8 = f8 = 64 bit float
    def __init__(self, data, seq_base=None, sf=None, precision='8'):
        self.data = data
        self.dtype = f"f{precision}"
        
        if seq_base is not None:
            if type(seq_base) is list:
                self.seq_base = seq_base
            else: 
                self.seq_base = BASES[seq_base]
        else:
            self.seq_base = list(sorted(set(data)))

        base_num = len(self.seq_base)
       
        if sf is not None:
            #todo:sf <= 1 sf >=1
            self.sf = sf
        else:
            self.sf =  1 - (sin(pi * (1/base_num)) / (sin(pi * (1/base_num)) + sin(pi * (1/base_num + 2*(floor(base_num/4) / base_num)))))
            
        if base_num == 4:
            corners = {
                "x":np.array([1,-1,-1,1], dtype=self.dtype),
                "y":np.array([-1,-1,1,1], dtype=self.dtype)
            }
        else:
            corners = self.__calc_corners()
        
        self.corners = pd.DataFrame.from_dict(corners, orient="index", columns=self.seq_base)
        
        self.coords = self.__calc_coords()
        
    def calc_fcgr(self, res=100):
        A = np.zeros((res,res), dtype=int)
        for i in range(len(self.data)):
            x = ceil((self.coords["x"][i]+self.r) * res/(2*self.r)) -1
            y = ceil((self.coords["y"][i]+self.r) * res/(2*self.r)) -1
            A[x,y] += 1
            
        return np.rot90(A)
    
    def __calc_coords(self):
        coords = {
            "x" : np.zeros(len(self.data), dtype=self.dtype),
            "y" : np.zeros(len(self.data), dtype=self.dtype)
        }
        last_x = 0
        last_y = 0
        for i, base in enumerate(self.data):
            corner = self.corners[base]
            last_x = last_x + (corner.x - last_x) * self.sf
            last_y = last_y + (corner.y - last_y) * self.sf
            coords["x"][i] = last_x
            coords["y"][i] = last_y
         
        return coords 
         
    def __calc_corners(self):
        base_num = len(self.seq_base)
        corners = {
            "x" : np.zeros(base_num, dtype=self.dtype),
            "y" : np.zeros(base_num, dtype=self.dtype)
        }
        for i in range(base_num):
            tmp = 2*i+1
            corners["x"][i] = self.r*sin(pi * (tmp/base_num))
            corners["y"][i] = self.r*cos(pi * (tmp/base_num))
        
        return corners

class RNAIntels():

    def __init__(self, args):
        self._fasta = args.filepath
        self._threads = args.threads

        self._model = load_model('RNAIntels.h5/RNAIntel.h5')
        data = self._read_fasta(self._fasta)
        identifiers, coded_data = self._encode_data(data) # multi-threading
        del data
        results = self._manage_prediction(identifiers, coded_data) # multi-threading

    def _read_fasta(self, filepath: str, min: int = 30, max: int = 60) -> object:
        sequences: list = []
        lengths: list = []
        name: list = []
        coding: object = pd.DataFrame()

        with open(filepath) as fasta_file:

            for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                sequences.append(str(seq_record.seq))
                lengths.append(len(seq_record.seq))
                name.append(seq_record.id)
                
            coding = pd.DataFrame({'sequence':sequences,'len':lengths,'id':name})

        # drop off sequences which are too long and to short
        coding_df = coding[coding['len'].between(min, max)]

        return coding_df

    def _encode_data(self, data: pd.DataFrame):
        # kick-out unknown nucleotids
        data['sequence'] = data['sequence'].map(lambda x: x.replace('N', ''))
        sequences = data.iloc[:, 0].values



        # multi-threading location
        coded_data: list = []
        identifiers: list = []
        counter: int = 0
        for x in sequences:
            cgr = CGR(data=x, seq_base=["A","G","T","C"])
            fcgr_matrix = cgr.calc_fcgr(res=8)
            coded_data.append(fcgr_matrix)
            identifiers.append(data['id'][counter])
            counter += 1

        return identifiers, coded_data

    def _manage_prediction(self, identifiers: list, coded_data: list):

        chunk_size: int = ceil(len(coded_data)/self._threads)
        seq_chunks: list = list(self.chunks(coded_data, chunk_size))
        id_chunks: list = list(self.chunks(identifiers, chunk_size))

        print(len(seq_chunks))

        #self._classify_sequences(id_chunks[0], seq_chunks[0])

        collection = Manager().list()
        sema = Semaphore(self._threads)
        
        jobs = []
        for i in range(0, len(seq_chunks)):
            p = Process(target=self._classify_sequences, args=(id_chunks[i], seq_chunks[i], sema))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()

        return collection


    def _classify_sequences(self, id_chunk: list, seq_chunk: list, sema):
        sema.acquire()

        if len(id_chunk) != len(seq_chunk):
            print("Not equal numbers of sequences and identifiers. Aborting")
            sema.release()
            exit(3)

        print("started chunks")

        print(len(id_chunk))
        print(len(seq_chunk))

        pred_probability = self._model.predict(np.array(seq_chunk))

        print(pred_probability)

        predictions = np.squeeze((pred_probability > 0.5).astype(int))

        print(predictions)

        results: list = []
        for iterator in range(0, len(predictions)):
            results.append((predictions[iterator], id_chunk[iterator]))

        print(results)

        sema.release()
        return True


    @staticmethod
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]


def main():
    parser = argparse.ArgumentParser(
                    prog='RNAIntels',
                    description='Classify RNA sequences into coding and noncoding',
                    epilog='Please cite us')
    parser.add_argument('filepath', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-t', '--threads', type=int, default=1)
    args = parser.parse_args()
    RNAIntels(args)

if __name__ == "__main__":
    main()
