#!/bin/sh
# used to;
# 1) create coreg isocontour image
# 2) compute the cost function of FLIRT coregistration
# 
# coreg_QC.sh <COREGISTERED_IMG> <REF(T1)_BRAIN> <REF(T1)_BRAIN_MASK> <OUTPUT_BASENAME>
# creates <OUTPUT_BASENAME>.png with isocontours image and <OUTPUT_BASENAME>.txt
# with costfunction value

in=$1
ref=$2
ref_mask=$3
out_base=$4

tmp=coreg_qc_tmp
mkdir -p $tmp

fslreorient2std $in $tmp/in.nii.gz
fslreorient2std $ref $tmp/ref.nii.gz
fslreorient2std $ref_mask $tmp/ref_mask.nii.gz

# sometimes this command simply produce a copy if nothing has to be changed.
# this causes a problem down the line if the file is nii instead of nii.gz...
num_unzipped_files=$(ls $tmp/*nii 2> /dev/null | wc -l)
if [ "$num_unzipped_files" != "0" ]; then
    unzipped_files=$(ls $tmp/*nii)
    for f in "${unzipped_files[@]}"
    do
        gzip $f
    done;
fi

cd $tmp
fslmaths in.nii.gz -nan in.nii.gz
fslmaths ref.nii.gz -nan ref.nii.gz
fslmaths ref_mask.nii.gz -nan ref_mask.nii.gz

# determine the lower and upper bounds for axial view
mkdir ax_slices && cd ax_slices
fslslice ../in.nii.gz img
i=1
ax_min=0
ax_max=0
for slice in `ls`
do
limits=$(fslstats $slice -R)
max=$(echo $limits | awk '{print $2}')
if [ "$max" != "0.000000" ] && [ $ax_min == 0 ]; then ax_min=$(($i+4)); fi
if [ "$max" == "0.000000" ] && [ $ax_min != 0 ]; then ax_max=$(($i-4)); break; fi
let i++
done
if [ $ax_max == 0 ]; then ax_max=$((i-1)); fi
cd .. && rm -rf ax_slices

# determine the sample space: we want roughly 7X5 axial slices
ax_range=$(($ax_max-$ax_min))
s=$(($ax_range/35))
ax_slices=`seq -$ax_max $s -$ax_min`

i=0
command=`echo pngappend ax_0.png`
for slice in $ax_slices
do
	slicer in.nii.gz ref.nii.gz -L -z $slice ax_$i.png
	let "z=$i%7"
	if [ "$i" != 0 ]; then
		if [ "$z" = 0 ]; then #linebreak
			command=`echo $command - ax_$i.png`
		else
			command=`echo $command + ax_$i.png`
		fi
	fi
	let i++
done
command=`echo $command ${out_base}.png`
$command

# Costfunction
flirt -in in.nii.gz -ref ref.nii.gz -refweight ref_mask.nii.gz -schedule $FSLDIR'/etc/flirtsch/measurecost1.sch' | head -1 | cut -f1 -d' ' > ${out_base}.txt
cd ../
cp $tmp/${out_base}.png ${out_base}.png
cp $tmp/${out_base}.txt ${out_base}.txt
rm -rf $tmp
