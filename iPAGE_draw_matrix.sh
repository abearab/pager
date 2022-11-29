function draw_matrix()
{
    export PAGEDIR='/data_gilbert/home/aarab/Workflows/iPAGE'
    
    mkdir -p ${1}_PAGE;

    local expfile=$1;
    local pvaluematrixfile=$2;
    local output=$3;
    
    perl ${PAGEDIR}/SCRIPTS/mi_go_draw_matrix.pl  \
            --pvaluematrixfile=${pvaluematrixfile} \
            --expfile=${expfile} \
            --order=1 --draw_sample_heatmap=false \
            --min=-3 --max=3 \
            --cluster=5 --quantized=0;
    wait

    mv -v ${expfile}_PAGE/${expfile}.summary.pdf $output
    rm -rv ${expfile}_PAGE/

}

draw_matrix $1 $2 $3
