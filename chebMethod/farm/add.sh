tag=$1
cd histos
types=`ls -1 *.root | grep -o "^[^_]*" | sort -u`
cd -
for t in $types
do
    hadd results/${tag}_${t}.root histos/${t}_*.root
done
