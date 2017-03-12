

import sys

with open(sys.argv[1]) as f:
    a = f.readlines()

with open(sys.argv[2]) as f:
    b = f.readlines()


def pdict(aa):
    aa = map(lambda x: x.strip().split(), aa)
    aadict = {}
    for x, y in aa:
        aadict[x] = float(y)
    return aadict

a = pdict(a)
b = pdict(b)

common = set(a.keys()) & set(b.keys())

for key in common:
    print "%s\t%.3f\t%.3f\t%.3f" % (key, a[key], b[key], a[key] / b[key])



