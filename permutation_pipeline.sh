



isoemfreq=$1
ec=$2
wdir=$3


mkdir -p $wdir
mkdir -p $wdir/permutations


for i in `seq 1 200`; do 
    mkdir -p $wdir/permutations/$i
    mkdir -p $wdir/permutations/$i/experiments
    mkdir -p $wdir/permutations/$i/results
    python permutation.py $isoemfreq 200 $wdir/permutations/$i/experiments
done



<<END
for exp in `ls $wdir/permutations/experiments`; do
    echo "RUNNING $exp";
    python em_ec_pathways.py /data1/igorm/SUPER_IGOR_ISOEM/isoem-1.1.5/stewart_initial/data/ec2path $ec $wdir/permutations/experiments/$exp $wdir/permutations/results/$exp
done
END





