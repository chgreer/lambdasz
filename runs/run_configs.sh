#!/bin/bash

#dir
dir=$2 #"/home/greer/sza_analysis/greer/mf/runs/sims/greer6"
MATOPT="-nosplash -nodesktop -singleCompThread -nodisplay "
counter=0

cd /home/greer/sza_analysis

NJOB=$1 #12

find $dir/*.in -print0 | while read -d $'\0' fnin
do
	counter=$((counter+1))
	fn=$(echo $fnin | sed 's/\(.*\).../\1/')
	echo "Running $fn"
	fnout=$fn.out
  fnlog=$fn.log
  nice matlab $MATOPT -r "mcmc('$fnin'); exit" > $fnlog &
	if [[ $counter -eq $NJOB ]] ; then
     echo "Reached $NJOB running jobs...waiting to complete."
     wait
     counter=0
     echo "Finished, starting new batch."
  fi

done
