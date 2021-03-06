
# to import rois

for name in 17_B05 39_B05
do
    /Applications/workbench/bin_macosx64/wb_command \
    -cifti-label-to-roi Conte69.parcellations_VGD11b.32k_fs_LR.dlabel.nii \
    Conte69.parcellations_VGD11b.32k_fs_LR.${name}.dscalar.nii -name ${name}
    /Applications/workbench/bin_macosx64/wb_command -cifti-separate Conte69.parcellations_VGD11b.32k_fs_LR.${name}.dscalar.nii COLUMN \
     -metric CORTEX_LEFT Conte69.parcellations_VGD11b.32k_fs_LR.${name}.dlabel.L.gii \
    -metric CORTEX_RIGHT Conte69.parcellations_VGD11b.32k_fs_LR.${name}.dlabel.R.gii
done


# to convert metric output files to labels, using .txt label definition file
for hemi in L R
do
    for t in pos neg
    do
        for pe in polar ecc
        do
            /Applications/workbench/bin_macosx64/wb_command \
            -metric-label-import ${hemi}_${t}_${pe}.metric.gii \
            ${hemi}_${t}_${pe}.metriclabel.txt \
            ${hemi}_${t}_${pe}.label.gii -unlabeled-value -1
        done
    done
done


