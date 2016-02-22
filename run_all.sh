#!/bin/bash -x

sequences="04 00 01 02 03 05 06 07 08 09 10"
depths="150 200"
inlier_thrs="3 4 5"
ransac_iters="300 400"

for sequence in $sequences
do
    for depth in $depths
    do
	for inlier_thr in $inlier_thrs
	do
	    for ransac_iter in $ransac_iters
	    do
		output_file="${sequence}_${depth}_${inlier_thr}_$ransac_iter"
		matlab -r "run_sequence('04','depth_thr',$depth,'inlier_thr',$inlier_thr,'ransac_iter',$ransac_iter);" > $output_file 2>&1 &
		echo "run $output_file"

	    done
	done
    done
    for job in `jobs -p`
    do
	echo $job
	wait $job || let "FAIL+=1"
    done
done
