#!/bin/bash 

Variable=(Run Lx_R Ly_R  BC_R S_R J_perp_R J_par_R H_R D_R Beta_R NT_R NSW_R NB_R NCPU_R)
Name=( Run Lx_R Ly_R  BC_R S_R J_perp_R J_par_R H_R D_R Beta_R NT_R NSW_R NB_R NCPU_R )
#      0    1    2     3    4      5      6      7   8     9    10   11   12    13
{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then
        export B_R_dir=`echo ${Variable[9]} | sed s/"\.0"//`
        
        export Dir="Lx"${Variable[1]}"_Ly"${Variable[2]}"_Jperp"${Variable[4]}"_Jpar"${Variable[5]}"_h"${Variable[7]}"_D"${Variable[8]}"_beta"$B_R_dir
        echo $Dir
        if [ ! -e $Dir ]; then
            mkdir $Dir
            cd  $Dir
            cp  ../Start/* . 
            cd ..
        else
            cd $Dir
            cp ../Start/parameters .
            cp ../Start/job.sh .
            cd ..
        fi
        cd $Dir
        
        let i=1
        while [  $i -lt 13 ]; do
            sed s/${Name[$i]}/${Variable[$i]}/    parameters  > tmp
            mv tmp parameters
            let i=i+1 
        done

        $SSE_DIR/src/ana.out *

        cd ..
    fi
    read -a Variable
    echo ${Variable[@]}
    echo ${Variable[0]}
done
}<sims

