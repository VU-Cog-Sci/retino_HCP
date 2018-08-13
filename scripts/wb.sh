
# to import rois

for name in BA2 BA1 BA3b BA4p BA3a BA4a BA6 BA17_V1 BA45 BA44 hOc5_MT BA18_V2 V3A V3B LO1 LO2 PITd PITv OP1 OP2 OP3 OP4 IPS1 IPS2 IPS3 IPS4 V7 V4v V3d 14c 13a 47s 14r 13m 13l 32pl 25 47m 47l Iai 10r 11m 11l 47r 10m Iam Ial 24 Iapm 10p V6_PHG06 ER 8_B05 6_B05 4_B05 9_B05 3_B05 1_B05 5_B05 7_B05 2_B05 31_B05 40_B05 44_B05 45_B05 23_B05 39_B05 43_B05 19_B05 47_B05 41_B05 30_B05 22_B05 42_B05 21_B05 38_B05 37_B05 20_B05 32_B05 24_B05 10_B05 25_B05 11_B05 46_B05 17_B05 18_B05 27_B05 36_B05 35_B05 28_B05 29_B05 26_B05 33_B05
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
