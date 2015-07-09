#!/bin/bash

#dir
dir="/home/greer/sza_analysis/greer/mf/runs/sims/greer9"
MATOPT="-nosplash -nodesktop -singleCompThread -nodisplay "
cd /home/greer/sza_analysis

#different runs
#REAL="01"
#nreal=300
SCAL="rozo"
DELTA="R500"
RHO="free"
COSM="lcdm"
NRUN=300

area="7398"
obs="ysz"

#fixed params
nstep=17500
burnin=2
burnsteps=1200
zmin=0.2
zmax=0.3
ncl=28
nreal=250

#obs param init guesses
init="1.73 1.58 0.15, 0.87, 5.3, 0.27, 0.0"

#rozo priors
priors="1.7, 0.87, 0.15, nan, nan, nan, nan"
prpiv=4.4
pr_wid="0.08, 0.18, 0.02, nan, nan, nan, nan"

for ((i=1; i<=$NRUN; i=i+1))
do
	if [ "$i" -le 9 ]
  then
    real=00$i
  elif [ "$i" -le 99 ]                                                          
  then
    real=0$i
  else
    real=$i
  fi
	for delta in $DELTA
	do
  	if [ "$delta" = "1Mpc" ]
  	then
			d=1
  	elif [ "$delta" = "R500" ]
  	then
			d=500
  	else
			d=2500
  	fi
  	for cosm in $COSM
  	do
    	for rho in $RHO
    	do
				if [ "$rho" = "fixed" ]
        	then
          	r=0
        	else
          	r=1
      	fi
				for scal in $SCAL
				do
					data_fn=/home/greer/sza_analysis/greer/maxbcg/maxbcg_ys_${scal}.txt
       		fn=sim_${cosm}_${delta}_${rho}rho_${scal}_tightprior_${real}real
	
					echo "Writing $fn.in"
					fn=${dir}/${fn}
	
      		#need code here to find pivot point 
					pivdir=/home/greer/sza_analysis/greer/mf/pivots
					pivot=${cosm}_${ncl}ncl_${zmin}-${zmax}zrange_${area}deg.pivot
      		if [ ! -e "$pivdir/$pivot" ] 
      		then
         		nice matlab $MATOPT -r "find_pivot_point($nreal,$ncl,$zmin,$zmax,$area,'$cosm'); exit" 
      		fi
      		pivot=$(cat $pivdir/$pivot)
					#pivot=6
	
      		fnin=$fn.in
      		fnout=$fn.out
      		echo "[options]" > $fnin
      		echo "nstep = $nstep" >> $fnin
      		echo "burnin = $burnin" >> $fnin
      		echo "burnsteps = $burnsteps" >> $fnin
      		echo "fn = $fnout" >> $fnin
      		echo "" >> $fnin
	
      		echo "[params]" >> $fnin
      		echo "init = $init" >> $fnin
      		echo "k = 0.70" >> $fnin
      		echo "fit_par = 1, 1, 1, 1, 1, 1, $r" >> $fnin
					echo "prior_cent_val = $priors" >> $fnin
					echo "prior_widths = $pr_wid" >> $fnin
      		echo "prior_pivot = $prpiv, $prpiv" >> $fnin
      		echo "mass_pivot = $pivot" >> $fnin
      		echo "flat_var_priors = 0,0" >> $fnin
      		echo "cosm = $cosm" >> $fnin
      		echo "" >> $fnin
	
      		echo "[obs]" >> $fnin
      		echo "type = $obs" >> $fnin
      		echo "zmin = $zmin" >> $fnin
      		echo "zmax = $zmax" >> $fnin
      		echo "area = $area" >> $fnin
      		echo "ncl = $ncl" >> $fnin
					echo "missing_y = 8,13,14,21,28" >> $fnin
					echo "rad = $d" >> $fnin
      		echo "issim = 1" >> $fnin
					#echo "data_fn = $data_fn" >> $fnin
    		done
			done
  	done
	done
done
