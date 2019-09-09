#!/bin/bash
# Modified version of create_spm_QC_images.sh to make it more generic to any 
# Background image (T1 or T2, or MNI space image). It requires GM and WM image
# in the same space as the background.
#
# In case the automatic bg_img intensity does not work well, it's possible to set 
# max val for BG image by specifying intensity percentile (<bg_max>). GM and WM 
# tissue maps will be thresholded at 0.001 and 1, but the max value can be set 
# by <gm_max> and <wm_max> (in this case bg_max needs to be inputted as well)
#
# Default is heat colormap, but can set use_jet to 1 to use jet colormap instead.
#
# create_vbm_QC_images.sh <bg_img> <gm_img> <wm_img> <use_jet> <out_name> <bg_max> <gm_max> <wm_max>

BG=$1
GM=$2
WM=$3
USE_JET=$4
OUTNAME=$5
BG_MAX=$6
GM_MAX=$7
WM_MAX=$8

# CREATION of QC images PNGs

export LANG=en_US
# this one is important, it sets the coding language in us
# otherwise seq will provide decimal numbers with a coma "," and not a point "."

tmp=vbm_qc_tmp
mkdir -p $tmp

fslreorient2std $BG $tmp/background.nii.gz
fslreorient2std $GM $tmp/gm_map.nii.gz
fslreorient2std $WM $tmp/wm_map.nii.gz

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

fslmaths background.nii.gz -nan background.nii.gz
fslmaths gm_map.nii.gz -nan gm_map.nii.gz
fslmaths wm_map.nii.gz -nan wm_map.nii.gz

# set overlay range if bg_min or bg_max is provided
if [ -z ${BG_MAX} ]; then 
    bg_range="-a";
else
    max_val=`fslstats background.nii.gz -P $BG_MAX`
    bg_range="0 ${max_val}";
fi

if [ -z ${GM_MAX} ]; then 
    gm_range="0.001 1";
else
    gm_range="0.001 ${GM_MAX}";
fi

if [ -z ${WM_MAX} ]; then 
    wm_range="0.001 1";
else
    wm_range="0.001 ${WM_MAX}";
fi

overlay_gm=`overlay 0 1 background.nii.gz $bg_range gm_map.nii.gz $gm_range gm_overlay.nii.gz`
echo $overlay_gm

overlay_wm=`overlay 0 1 background.nii.gz $bg_range wm_map.nii.gz $wm_range wm_overlay.nii.gz`
echo $overlay_wm

# determine the lower and upper bounds for native gm/wm maps for both 3 views

## for axial view
mkdir ax_slices && cd ax_slices
fslslice ../gm_map.nii.gz gm
i=1
ax_min=0
ax_max=0
for slice in `ls`
do
limits=$(fslstats $slice -R)
max=$(echo $limits | awk '{print $2}')
var=$(awk 'BEGIN{ print "0.1"<"'$max'" }')
if [ "$var" -eq 1 ] && [ $ax_min == 0 ]; then ax_min=$(($i)); fi
if [ "$var" -ne 1 ] && [ $ax_min != 0 ]; then ax_max=$(($i)); break; fi
let i++
done
cd .. && rm -rf ax_slices

# Find interval that will show at least 30 slices
ax_interval=$((($ax_max - $ax_min)/30))
ax_start=$(($ax_max-$ax_interval))
ax_end=$(($ax_min+$ax_interval))
ax_slices=`seq -$ax_start $ax_interval -$ax_end`

## for coronal view
# Get the x and y dim size
gm_info=$(fslinfo gm_map.nii.gz)
xdim=$(echo $gm_info | awk '{ print $4}')
ydim=$(echo $gm_info | awk '{ print $6}')

fslswapdim gm_map.nii.gz x z -y cor.nii.gz
mkdir cor_slices && cd cor_slices
fslslice ../cor.nii.gz gm
i=1
cor_min=0
cor_max=0
for slice in `ls`
do
limits=$(fslstats $slice -R)
max=$(echo $limits | awk '{print $2}')
var=$(awk 'BEGIN{ print "0.1"<"'$max'" }')
if [ "$var" -eq 1 ] && [ $cor_min == 0 ]; then cor_min=$(($i)); fi
if [ "$var" -ne 1 ] && [ $cor_min != 0 ]; then cor_max=$(($i)); break; fi
let i++
done

# Take 7x5 slices to display
cor_interval=$((($cor_max - $cor_min)/30))

## rearrangment to match with the swapping
new_cor_min=$(($ydim-$cor_max+$cor_interval))
new_cor_max=$(($ydim-$cor_min-$cor_interval))
cd .. && rm -rf cor_slices
rm cor.nii.gz

cor_slices=`seq -$new_cor_max $cor_interval -$new_cor_min`

## for sagittal view
fslswapdim gm_map.nii.gz z y -x sag.nii.gz
mkdir sag_slices && cd sag_slices
fslslice ../sag.nii.gz gm
i=1
sag_min=0
sag_max=0
for slice in `ls`
do
limits=$(fslstats $slice -R)
max=$(echo $limits | awk '{print $2}')
var=$(awk 'BEGIN{ print "0.1"<"'$max'" }')
if [ "$var" -eq 1 ] && [ $sag_min == 0 ]; then sag_min=$(($i)); fi
if [ "$var" -ne 1 ] && [ $sag_min != 0 ]; then sag_max=$(($i)); break; fi
let i++
done

sag_interval=$((($sag_max - $sag_min)/30))
## rearrangment to match with the swapping
new_sag_min=$(($xdim-$sag_max+$sag_interval))
new_sag_max=$(($xdim-$sag_min-$sag_interval))
cd .. && rm -rf sag_slices
rm sag.nii.gz

sag_slices=`seq -$new_sag_max $sag_interval -$new_sag_min`


maps="gm_overlay.nii.gz wm_overlay.nii.gz"

if [ "$USE_JET" == 0 ]; then
    color_lut="";
else
    color_lut="-l ${FSLDIR}/etc/luts/renderjet.lut ";
fi

for map in $maps
do
	name=$(echo $map | cut -d '.' -f1)
    
	# sagittal
	i=0
	command=`echo pngappend sag_0.png`
	
	for slice in $sag_slices
		do
			slicer $map -L ${color_lut}-x $slice sag_$i.png
			let "z=$i%7"
			if [ "$i" != 0 ]; then
				if [ "$z" = 0 ]; then #linebreak
					command=`echo $command - sag_$i.png`
				else
					command=`echo $command + sag_$i.png`
				fi
			fi
			let i++
		done
	
	command=`echo $command ${OUTNAME}_${name}_sag.png`
	$command
	rm sag_*

	# coronal
	i=0
	command=`echo pngappend cor_0.png`

	for slice in $cor_slices
		do
			slicer $map -L ${color_lut}-y $slice cor_$i.png
			let "z=$i%7"
			if [ "$i" != 0 ]; then
				if [ "$z" = 0 ]; then #linebreak
					command=`echo $command - cor_$i.png`
				else
					command=`echo $command + cor_$i.png`
				fi
			fi
			let i++
		done
	
	command=`echo $command ${OUTNAME}_${name}_cor.png`
	$command
	rm cor_*

	# axial
	i=0
	command=`echo pngappend ax_0.png`
	
	for slice in $ax_slices
		do
			slicer $map -L ${color_lut}-z $slice ax_$i.png
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
	
	command=`echo $command ${OUTNAME}_${name}_ax.png`
	$command
	rm ax_*
done

cd ../
cp $tmp/${OUTNAME}*png .
rm -rf $tmp
