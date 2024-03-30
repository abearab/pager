function draw_matrix()
{
    if [ -z ${TEISERDIR+x} ]; then
        echo "Error: TEISERDIR is not set. Please set TEISERDIR to the TEISER directory.";
        return;
    fi
    
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
