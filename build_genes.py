

import sys
from Bio import SeqIO
import re



fasta = sys.argv[1]
gtf = sys.argv[2]
output = sys.argv[3]




sw = {}
handle = open(fasta, "rU")
for record in SeqIO.parse(handle, "fasta"):
    sw[record.id] = record.seq.tostring()
handle.close()



with open(gtf) as f:
    a = f.readlines()

a = map(lambda x: x.strip().split("\t"), a)



def reverse_complement(seq):    
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases


def either_revcomp(seq, strand):
    if strand == "-1":
        return reverse_complement(seq)
    else:
        return seq


resulting_fasta_genes = {}
for line in a:
    chromosome_gtf = line[0]
    tag = re.search("locus_tag=(.+);", line[-1])
    if tag:
        locus_tag = tag.group(1)
        start, end, strand = int(line[3]) - 1, int(line[4]) - 1, line[6]
        resulting_fasta_genes[locus_tag] = either_revcomp(sw[chromosome_gtf][start: end + 1], strand)



with open(output, "w") as f:
    for x, y in resulting_fasta_genes.items():
        f.write(">%s\n" % x)
        f.write("%s\n" % y)



