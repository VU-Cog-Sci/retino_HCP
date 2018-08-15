# 1. open freeview

export sub_path=/Users/martin/disks/ae_S/2018/visual/nprf_hcp/pp/freesurfer/fsaverage
cd $sub_path
freeview -v mri/T1.mgz \
mri/wm.mgz:visible=0 \
mri/brainmask.mgz:colormap=heat:opacity=0.5 \
-f surf/lh.white:edgecolor=yellow \
surf/rh.white:edgecolor=yellow \
surf/lh.pial:edgecolor=red \
surf/rh.pial:edgecolor=blue \
surf/rh.inflated:visible=0 \
surf/lh.inflated:visible=0 

# in aeneas (first set your SUBJECTS_DIR in bash profile to where are your subject to recon)
tmux
recon-all -autorecon-pial -subjid sub-002
#recon-all -autorecon-pial -subjid sub-003
#recon-all -autorecon-pial -subjid sub-004


# 2. cut the brain for flattening it (https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOccipitalFlattenedPatch)
# template http://gallantlab.org/huth2016/

# 1. Start cutting at the calcarine gyrus (at the edge of V1) and follow it until the first sulcus appears, then go straight until the midbrain
# Click on show atlas
# 2. Start cutting between precuneus and paracentral and parietal superior. Go from that point along the precuneus down to the center
# 3. Draw line from corner of lateral view (pink), passing between cingular medial posterior (blue) and cingular medial anterior (dark blue)
# 4. start at intersection front superior gyrus / front middle gyrus, draw line passing through trans_frontopol., go along the line of the recturs down to the center of the brain
# 5. Start at edge of temporal pole, passing throug temporal pole(yellow), temporal superior plan polar (light green), then pass in the middle of the temporal inferior (light_yellow)
# 6. cut from the top the corpus callosum including the subcalosum and the non anottaded area (including hypocampal)

# on laptop
tksurfer sub-003 lh inflated -curv -annotation aparc.a2009s
tksurfer sub-003 rh inflated -curv -annotation aparc.a2009s

# >> save as lh.full.patch.3d 
# >> save as rh.full.patch.3d 

# run on aeaneas the flattening function
# cd to surface folder
tmux
cd /home/shared/2018/visual/nprf_hcp/pp/freesurfer/fsaverage/surf/
mris_flatten lh.full.patch.3d lh.full.flat.patch.3d
mris_flatten rh.full.patch.3d rh.full.flat.patch.3d

# to check
tksurfer sub-004 lh inflated -patch lh.full.flat.patch.3d -curv