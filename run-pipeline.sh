#!/bin/bash


reads1=$1
reads2=$2
Index=$3
GTF=$4
fna=$5
outdir=$6
isoem2=$7



hisat2-build -f $fna $Index


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
       -O $outdir/isoem2_results \
       -C 95 \
       $outdir/hisat2.sam ) \
       2> $outdir/isoemlog \
       > $outdir/log



