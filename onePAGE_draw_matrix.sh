function draw_matrix()
{
    export TEISERDIR='/data_gilbert/home/aarab/Workflows/TEISERv1.1'
    
    local pdir=$1;
    local exp=$2;
    
    wd = `pwd`
    cd $pdir
    perl ${TEISERDIR}/Scripts/teiser_draw_matrix.pl \
        --pvmatrixfile=${exp}.txt.matrix \
        --summaryfile=${exp}.txt.summary \
        --expfile=${exp}.txt \
        --quantized=0 \
        --order=0 \
        --min=-3 --max=3 --cluster=5
    wait    
    
    cd $wd
}

draw_matrix $1 $2
