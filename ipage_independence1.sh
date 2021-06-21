## $1: input file contain table with first two columns as gene name/id and numeric value
species='human'
# pattern='human_*gs*'
pattern='human_ensembl*'

expfile=`basename $1`
outdir=${1/.txt/};
ipage_ann='/flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS/'


mkdir -p $outdir
cd $outdir; cd ../

echo `pwd`

for f in `ls -d ${ipage_ann}/${pattern}`; do

    base=`basename "$f"`;
    echo '________________' $base '________________';

    if [ -d "${outdir}/${base}/" ]; 
    then
        echo 'This result exist!';
    ## TO-DO; Force run to remove stuff from previous run. 
    else
        # Run iPAGE 
        perl $PAGEDIR/page.pl --expfile=$expfile --independence=1 \
        --species=$base --exptype=continuous --ebins=11 --nodups=1; 
        wait

        mv -v ${expfile}_PAGE/ ${outdir}/${base}/;

        # remove the result folder if it was empty
        counter="$(wc -l < ${outdir}/${base}/pvmatrix.txt)"
        if [ $counter -le "$(echo '1')" ] || [ -z $counter ]
        then
            rm -r ${outdir}/${base}/
        else
            # keep complete results 
            pv=${outdir}/${base}/pvmatrix.txt;
            pv0=${pv/.txt/.all.txt}; pvL=${pv/.txt/.L.txt}; pvR=${pv/.txt/.R.txt};

            # subset pathways with p-value > 2 in the first (Left) or last (Right) cluster
            mv -v $pv $pv0;
            cat $pv0 | awk -F'\t' 'BEGIN{FS=OFS="\t"}; $2>2 {print $0}' > $pvL;
            cat $pv0 | awk -F'\t' 'BEGIN{FS=OFS="\t"}; $12>2 {print $0}' > $pvR;

            for sum in `ls ${outdir}/${base}/*.summary*`; do 
                sum0=${sum/summary/all};
                mv -v $sum $sum0;
            done
            wait

            # draw sided heatmaps     
            declare -a Sides=('L' 'R');

            for side in "${Sides[@]}"; do

                # include if pv matrix contain any values 
                counter="$(wc -l < ${outdir}/${base}/pvmatrix.${side}.txt)"
                if [ $counter -le "$(echo '1')" ] || [ -z $counter ]
                then 
                    echo '_____...________' $base $side '-> No signal!';    
                else 
                    echo '_____...________' $base $side '...';
                    mkdir -p ${expfile}_PAGE;
                    perl $PAGEDIR/SCRIPTS/mi_go_draw_matrix.pl  \
                    --pvaluematrixfile=${outdir}/${base}/pvmatrix.${side}.txt \
                    --expfile=$expfile \
                    --order=1 --draw_sample_heatmap=false \
                    --min=-3 --max=3 \
                    --cluster=5 --quantized=0;
                    wait
                    for sum in `ls ${expfile}_PAGE`; do 
                        sumS=${sum/summary/$side};
                        mv -v ${expfile}_PAGE/${sum} ${outdir}/${base}/${sumS};
                    done     
                    rm -vr ${expfile}_PAGE;
                fi;

            done
            
            mv -v $pv0 $pv
            
            # pdf2png!
            for file in ${outdir}/${base}/${expfile}*.pdf; do 
                file=${file/.pdf/}
                bash /rumi/shams/abe/Workflows/my_scripts/pdf2png.sh ${file}.pdf
                new=${file/\/$expfile/}
                cp -v ${file}.pdf ${new}.pdf
                cp -v ${file}.png ${new}.png
            done 
            
        fi;
        
    fi

done
