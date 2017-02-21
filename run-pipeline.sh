#!/bin/bash


reads1=$1
reads2=$2
Index=$3
GTF=$4
fna=$5
outdir=$6
isoem2=$7



#isoem2=../../../../../bin/isoem2 #/data1/sahar/isoem2experiments/pipeline2/code/isoem-b/bin/isoem2
#Index=index/index
#GTF=89935.assembled.fna.GTF

#hisat2-build -f 89935.assembled.fna index/index
hisat2-build -f $fna $index


(time \
   /usr/local/bin/hisat2 \
        -p 30 \
        --reorder --no-discordant --no-mixed --no-unal --no-hd --sensitive \
        -x ${Index} \
        -1 ${reads1} -2 ${reads2} ) \
        2> $outdir/hisat2log> $outdir/hisat2.sam


(time ${isoem2} \
       -a \
       -G ${GTF} \
       -O $oudir/isoem2_results \
       .$outdir/hisat2.sam ) \
       2> $outdir/isoemlog \
       > $outdir/log



