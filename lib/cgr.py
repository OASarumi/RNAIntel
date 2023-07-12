import numpy as np
import pandas as pd

from math import ceil

# This code is based on parts of the R package kaos by Dominic Eger and Hannah Franziska Löchel, implemented in Python by Hannah Franziska Löchel and Sandra Clemens
class CGR:
    r = 1

    #preciion 8 = f8 = 64 bit float
    def __init__(self, data, precision='8'):
        self.data = data
        self.dtype = f"f{precision}"
        self.seq_base=["G","T","C","A"]
        self.sf=0.5
        base_num = 4


        corners = {
            "x":np.array([1,-1,-1,1], dtype=self.dtype),
            "y":np.array([-1,-1,1,1], dtype=self.dtype)
        }

        self.corners = pd.DataFrame.from_dict(corners, orient="index", columns=self.seq_base)


        self.coords = self.__calc_coords()

    def calc_fcgr(self, res=100):
        A = np.zeros((res,res), dtype=int)
        for i in range(len(self.data)):
            x = ceil((self.coords["x"][i]+self.r) * res/(2*self.r)) -1
            y = ceil((self.coords["y"][i]+self.r) * res/(2*self.r)) -1
            A[x,y] += 1

        return A

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
