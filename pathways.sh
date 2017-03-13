


isoem2=/data1/igorm/SUPER_IGOR_ISOEM/isoem-1.1.5/bin/isoem2
hisat2=/usr/local/bin/hisat2
fastaToGTF=`dirname $isoem2`/fastaToGTF
ec2path=/data1/igorm/SUPER_IGOR_ISOEM/isoem-1.1.5/stewart_initial/data/ec2path
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


wdir=$1
reads=$2 # must be fastq
contigs=$3 # must be fasta
gff=$4 # annotation with genes
ec=$5 # these are ec numbers mapped



mkdir -p $wdir



echo "SPLITTING READS"
python $DIR/scripts/split.fastq.py $reads $wdir
echo "DONE"


python $DIR/scripts/build_genes.py $contigs $gff $wdir/metacontigs.fasta

$fastaToGTF $wdir/metacontigs.fasta


reads1=$wdir/r1.fastq
reads2=$wdir/r2.fastq
Index=$wdir/hisat2_index
GTF=$wdir/metacontigs.fasta.GTF
fna=$wdir/metacontigs.fasta
outdir=$wdir/results


mkdir -p $Index
mkdir -p $outdir

${hisat2}-build -f $fna $Index


(time \
   $hisat2 \
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






python $DIR/scripts/em_ec_pathways.py $ec2path $ec $wdir/results/isoem2_results/hisat2/output/Genes/gene_fpkm_estimates $wdir/metabolic.pthws



mkdir -p $wdir/bootstrap


for i in `seq 1 200`; do
    python $DIR/scripts/em_ec_pathways.py $ec2path $ec $wdir/results/isoem2_results/hisat2/bootstrap/experiment_${i}/Genes/gene_fpkm_estimates $wdir/bootstrap/metabolic.pthws.${i}
done






