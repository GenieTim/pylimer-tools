# Using the generator pattern (an iterable)
from pylimer_tools.io.readLammpData import readLammpData

"""
This class represents a sequence of Universes, 
with the Universe's data files only being read on access.
"""


class LazyUniverseDataFileSequence(object):

    def __init__(self, files):
        self.files = files
        self.n = len(self.files)
        self.idx = 0
        self.parsedFiles = {}

    def parseFile(self, idx):
        self.parsedFiles[idx] = readLammpData(self.files[idx])
        return self.parsedFiles[idx]

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        if self.idx < self.n:
            cur, self.idx = self.idx, self.idx+1
            if (cur in self.parsedFiles):
                return self.parsedFiles[cur]
            else:
                return self.parseFile(cur)
        raise StopIteration()
