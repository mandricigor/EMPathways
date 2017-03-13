

import sys

tmpfile=sys.argv[1]
wdir = sys.argv[2]

fastq_1 = open('%s/r1.fastq' % wdir, "w")

fastq_2 = open('%s/r2.fastq' % wdir, "w")

[fastq_1.write(line) if (i % 8 < 4) else fastq_2.write(line) for i, line in enumerate(open('%s' % tmpfile))]

fastq_1.close()

fastq_2.close()
