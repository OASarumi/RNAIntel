import os
import numpy as np
import pandas as pd
import time
import gzip

from .cgr import CGR
from multiprocessing import Semaphore, Manager, Process
from Bio import SeqIO
from math import ceil


class RNAIntels():

    def __init__(self, args):
        self._seq_file = args.sequence_file
        self._threads = args.threads
        self._verbose = args.verbose
        self._model_path = args.model
        self._start_time = time.time()

        # read sequence file
        data = self._read(self._seq_file, args.min, args.max)

        # overwrite original data with matrices
        data = self._manage_encoding(data)

        # perform model classification
        results = self._classify_sequences(data)

        # format results output
        self._print_results(results)

        if self._verbose >= 1:
            print("# Finished (%d s)" % (time.time()-self._start_time))

    def _read(self, filepath: str, min: int = 30, max: int = 60) -> pd.DataFrame:
        sequences: list = []
        lengths: list = []
        name: list = []
        coding: object = pd.DataFrame()

        if self._verbose >= 1:
            print("# Reading sequence file")

        # Check whetever file is gzip compressed or not
        if filepath[-5:len(filepath)].lower() == ".gzip" or filepath[-3:len(filepath)].lower() == ".gz":
            if self._verbose >= 1:
                print("## GZIP compression assumed")
            handler = gzip.open(filepath, "rt")
            remaining_name = filepath[0:-5] if filepath[-5:len(filepath)] == ".gzip" else filepath[0:-3]
        else:
            handler = open(filepath, "r")
            remaining_name = filepath

        # Check for FASTA/FASTQ
        if remaining_name[-6:len(remaining_name)].lower() == ".fastq":
            if self._verbose >= 1:
                print("## Detect FASTQ file")
            sequence_type = "fastq"
        else:
            sequence_type = "fasta"

        with handler as sequence_file:

            for seq_record in SeqIO.parse(sequence_file, sequence_type):
                sequences.append(str(seq_record.seq))
                name.append(seq_record.id)
                
            data = pd.DataFrame({'sequence':sequences, 'id':name})

        # drop off sequences which are too long and to short
        between_filter = lambda x: len(x) >= min and len(x) <= max
        data: pd.DataFrame = data.loc[data['sequence'].apply(between_filter)]
        data['sequence'] = data['sequence'].map(lambda x: x.replace('N', ''))

        return data.reset_index(drop=True)

    @staticmethod
    def _get_split_ranges(total_len: int, divider: int):
        chunk_size: int = ceil(total_len/divider)
        counter: int = 0
        ranges: list = []
        while counter <= total_len:
            tmp: int = counter
            counter += chunk_size
            ranges.append([tmp, counter])

        if ranges[-1][1] > total_len:
            ranges[-1][1] = total_len

        return ranges

    def _manage_encoding(self, data: pd.DataFrame) -> list:
        # Multi-threaded mode
        if self._threads > 1:
            if self._verbose >= 1:
                print("# Encoding: Multi-threaded mode (%d s)" % (time.time()-self._start_time))

            # Resuts bucket and semaphore
            collection = Manager().list()
            sema = Semaphore(self._threads)

            # Indeces to split table into mostly equally long packages
            ranges = self._get_split_ranges(len(data), self._threads)

            jobs = []
            for i in range(0, len(ranges)):
                seqs = data["sequence"][range(ranges[i][0], ranges[i][1])].values
                ids = data["id"][range(ranges[i][0], ranges[i][1])].values
                p = Process(target=self._encoding_mt, args=(seqs, ids, sema, collection))
                jobs.append(p)
                p.start()
            for proc in jobs:
                proc.join()
            encoded_data = collection

        # Single-threaded mode
        else:
            if self._verbose >= 1:
                print("# Encoding: Single-threaded mode (%d s)" % (time.time()-self._start_time))
            encoded_data = self._encoding_st(data["sequence"], data["id"])

        return encoded_data

    def _encoding_st(self, seqs: list, ids: list):
        results: list = []

        for i in range(0, len(seqs)):
            try:
                fcgr = CGR(data=seqs[i], seq_base=["A","G","T","C"]).calc_fcgr(res=8)
                results.append((ids[i], fcgr))
            except:
                continue

        return results

    def _encoding_mt(self, seqs: list, ids: list, sema: Semaphore, collection):
        sema.acquire()

        for i in range(0, len(seqs)):
            try:
                fcgr = CGR(data=seqs[i], seq_base=["A","G","T","C"]).calc_fcgr(res=8)
                collection.append((ids[i], fcgr))
            except:
                continue

        sema.release()
        return True

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

    def _classify_sequences(self, data: list):
        from tensorflow.keras.models import load_model

        if self._model_path is not None:
            model = load_model(self._model_path)
        else:
            __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
            model = load_model(os.path.join(__location__, '../RNAIntelsModel.h5'))

        if self._verbose >= 1:
            print("# Classification (%d s)" % (time.time()-self._start_time))

        matrices = np.array([x[1] for x in data])
        pred_probability = model.predict(matrices)
        predictions = np.squeeze((pred_probability > 0.5).astype(int))

        if self._verbose >= 1:
            print("# positive predictions: %d / %d" % (sum(predictions), len(predictions)) )

        results: list = []
        for i in range(0, len(predictions)):
            results.append((data[i][0], predictions[i]))

        return results

    def _print_results(self, results: list):
        for result in results:
            print("%s\t%d" % (result[0], result[1]))


    @staticmethod
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]