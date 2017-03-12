




ec_pathways_dict=ec_pathways_dict.txt
isoem2=/data1/igorm/SUPER_IGOR_ISOEM/isoem-1.1.5/bin/isoem2
fastaToGTF=../../bin/fastaToGTF


wdir=$1



merged_fastq=`ls $wdir/$wdir/RawData/*.fastq | head -1`

python split.fastq.py $merged_fastq $wdir/$wdir/RawData/



gomion=`ls $wdir/GoMion_* | head -1`


fna=`ls $wdir/$gomion/IMGData/*.fna | head -1`

$fastaToGTF $fna



reads1=$wdir/$wdir/RawData/r1.fastq
reads2=$wdir/$wdir/RawData/r2.fastq
Index=$wdir/hisat2_index/index
mkdir -p $wdir/hisat2_index
GTF=${fna}.GTF
fna=$fna
outdir=$wdir/results
mkdir -p $outdir


bash run-pipeline.sh $reads1 $reads2 $Index $GTF $fna $outdir $isoem2


ec=`ls $wdir/$gomion/IMGData/*.EC | head -1`

python em_ec_pathways.py $ec_pathways_dict $ec $outdir/isoem2_results/hisat2/output/Genes/gene_fpkm_estimates $wdir/pathway_abundance.txt








