


import sys
from scipy.stats import pearsonr


file1 = sys.argv[1]
file2 = sys.argv[2]


def readf(filex):
    with open(filex) as f:
        a = f.readlines()
    a = map(lambda x: x.strip().split(), a)
    adict = {}
    for x, y in a:
        adict[x] = float(y)
    return adict



v1 = readf(file1)
v2 = readf(file2)


common = set(v1.keys()) & set(v2.keys())

a = []
b = []
for key in common:
    a.append(v1[key])
    b.append(v2[key])


print pearsonr(a, b)







