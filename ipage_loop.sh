## $1: input file contain table with first two columns as gene name/id and numeric value
pattern='msigdb*'
# pattern='human_*_gs*'

expfile=`basename $1`
outdir=${1/.txt/};

mkdir -p $outdir
cd $outdir; cd ../

for f in `ls -d $PAGEDIR/PAGE_DATA/ANNOTATIONS/${pattern}`; do
    base=`basename "$f"`;
    echo '________________' $base '________________';
    if [ -d "${outdir}/${base}/" ]; 
    then
        echo 'This result exist!';
    else
        # Run iPAGE 
        perl $PAGEDIR/page.pl --expfile=$expfile --species=$base --exptype=continuous --ebins=11 --nodups=1             --independence=0;
        # --independence=0; option for comparing results between multiple smaples.
        wait
        mv -v ${expfile}_PAGE/ ${outdir}/${base}/;
        # remove the result folder if it was empty
        counter="$(wc -l < ${outdir}/${base}/pvmatrix.txt)"
        if [ $counter -le "$(echo '1')" ] || [ -z $counter ]
        then
            rm -r ${outdir}/${base}
        fi;
        rm -r \
          ${outdir}/${base}/*txt.summary* \
          ${outdir}/${base}/SUMMARY \
          ${outdir}/${base}/imgsrc \
          ${outdir}/${base}/info.txt \
          ${outdir}/${base}/pvmatrix.txt.killed \
          ${outdir}/${base}/pvmatrix.txt.log
          # ${outdir}/${base}/input.ipage_quantized
    fi
done
