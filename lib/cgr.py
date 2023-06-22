import numpy as np
import pandas as pd

from math import sin, cos, pi, floor, ceil

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