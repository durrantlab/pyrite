import copy
from collections import OrderedDict

def array(lst):
    if len(lst) == 0:
        return Array1D(lst)
    elif type(lst[0]) is list:
        return Array2D(lst)
    else:
        return Array1D(lst)


class Array:
    lst = []
    shape = ()

    def __iter__(self):
        # So you could use for... in...
        for x in self.lst:
            yield x

    def __str__(self):
        return str(self.lst)

    def __repr__(self):
        return "array(" + self.__str__() + ")"

    def __getitem__(self, selection):
        # []
        new_lst = []

        try:
            # Assume the selection is iterable.
            for i in selection:
                new_lst.append(self.lst[i])
        except:
            # selection must not be iterable. Just a number.
            new_lst.append(self.lst[selection])

        # If it's just one item, then just return the value, not the value in
        # an array
        if len(new_lst) == 1:
            return new_lst[0]

        return array(new_lst)
    
    def __setitem__(self, key, item):
        try:
            # Assume the selection (key) is iterable.
            for i in key:
                self.lst[i] = item
        except:
            # selection must not be iterable. Just a number?
            self.lst[key] = item
    
    def __len__(self):
        return len(self.lst)


class Array1D(Array):
    def __init__(self, lst):
        self.shape = (len(lst),)
        self.lst = lst

    def __eq__(self, other):
        bools = copy.deepcopy(self.lst)
        for x in range(self.shape[0]):
            bools[x] = (bools[x] == other)
        return array(bools)
    

class Array2D(Array):
    def __init__(self, lst):
        self.shape = (len(lst), len(lst[0]))

        self.lst = []
        for row in lst:
            self.lst.append(array(row))

    def __eq__(self, other):
        bools = copy.deepcopy(self.lst)
        for x in range(self.shape[0]):
            for y in range(self.shape[1]):
                bools[x][y] = (bools[x][y] == other)
        return array(bools)
    

class DType:
    names = []
    descr = []


class DictArray:
    # Collection of arrays accessible through dictionary keys.
    dict = OrderedDict({})
    dtype = DType()
    ndim = 1  # I think for a DictArray this is always 1?

    def __init__(self, dict, dtypes = []):
        # dict maps keys to arrays.
        for key in dict:
            self.dict[key] = array(dict[key])
        
        self.dtype.names = dict.keys()
        for i, key in enumerate(dict.keys()):
            try:
                self.dtype.descr.append((key, "|" + dtypes[i]))
            except:
                pass
    
    def __repr__(Self):
        return "dictarray(" + self.__str__() + ")"
    
    def __str__(self):
        return str(self.dict)

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return self.dict[key]
        else:
            # So key must be something that should act on the individual
            # Arrays, like [1, 2, 3, 4]
            new_dict = OrderedDict({})
            for str_key in self.dict.keys():
                new_dict[str_key] = self.dict[str_key][key].lst
            
            return DictArray(new_dict)

