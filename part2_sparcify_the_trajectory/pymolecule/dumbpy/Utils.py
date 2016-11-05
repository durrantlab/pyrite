from Array import DictArray
from Array import array
from collections import OrderedDict

def genfromtxt(fname, dtype = "", names = [], delimiter = []):
    # fname here is a file object.

    # Load the data
    data = OrderedDict({})
    for n in names:
        data[n] = []

    for line in fname:
        parts = []
        for num in delimiter:
            parts.append(line[:num])
            line = line[num:]
        
        for i, name in enumerate(names):
            data[name].append(parts[i])

    # Fix any types
    dtype = [l.strip() for l in dtype.split(",")]
    for i, name in enumerate(names):
        if "int" in dtype[i] or dtype[i][:1] == "i":
            data[name] = [int(l) for l in data[name]]
        if "float" in dtype[i] or dtype[i][:1] == "f":
            data[name] = [float(l) for l in data[name]]
    
    return DictArray(data, dtype)

def nonzero(arr):
    if len(arr.shape) == 1:
        indx_to_keep = []
        for i, val in enumerate(arr):
            if val != 0:
                indx_to_keep.append(i)
        indx_to_keep = (array(indx_to_keep),)
    elif len(arr.shape) == 2:
        indx1_to_keep = []
        indx2_to_keep = []

        for i, row in enumerate(arr):
            for j, val in enumerate(row):
                if val != 0:
                    indx1_to_keep.append(i)
                    indx2_to_keep.append(j)
        indx_to_keep = (array(indx1_to_keep), array(indx2_to_keep))

    #print arr
    #print arr.shape
    #sdf
    return indx_to_keep

def logical_or(arr1, arr2):
    if len(arr1.shape) == 1:
        or_result = [x or y for x,y in zip(arr1, arr2)]
        return array(or_result)

def defchararray_strip(arr):
    if len(arr.shape) == 1:
        strip_result = [s.strip() for s in arr]
        return array(strip_result)
    