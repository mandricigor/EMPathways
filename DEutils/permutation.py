


import sys
import random

freqfile = sys.argv[1]
nperm = sys.argv[2]
outdir = sys.argv[3]


with open(freqfile) as f:
    a = f.readlines()

freq = map(lambda x: x.strip().split(), a)

features = map(lambda x: x[0], freq)
values = map(lambda x: float(x[1]), freq)




for i in range(1, int(nperm) + 1):
    random.shuffle(values)
    with open(outdir + "/" + "perm_%s" % i, "w") as f:
        for x, y in zip(features, values):
            f.write("%s\t%s\n" % (x, y))



