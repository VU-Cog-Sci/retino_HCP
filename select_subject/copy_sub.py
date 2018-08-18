"""
-----------------------------------------------------------------------------------------
copy_sub
-----------------------------------------------------------------------------------------
Goal of the script:
create jobscript to run in a cluster (LISA or CARTESIUS) or server
-----------------------------------------------------------------------------------------
Input(s):
none
-----------------------------------------------------------------------------------------
Output(s):
none
-----------------------------------------------------------------------------------------
Exemple:
# to run from lisa
cd /Users/martin/Dropbox/GitHub/retino_HCP/ 
python select_subject/copy_sub.py
-----------------------------------------------------------------------------------------
"""

import os
import glob
import ipdb

subs = ['192641','105923','111312','926862','182739','167440','789373','690152','148133','467351',\
		'233326','144226','638049','318637','165436','771354','995174','463040','187345','547046',\
		'150423','167036','104416','536647','397760','552241','169040','973770','116726','130114']

files_to_trans = [	'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries_L.func_bla_psc_av.gii',\
					'tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries_R.func_bla_psc_av.gii'	]

database_dir = '/home/shared/2018/visual/nprf_hcp/data/'
target_dir = 'szinte@lisa.surfsara.nl:/home/szinte/data/nprf_hcp'

copy_cmd_beg = 'rsync -a --no-g --no-p -vzhe ssh --progress'

for sub in subs:
	for file_to_trans in files_to_trans:
		from_file = '%s/%s/%s'%(database_dir,sub,file_to_trans)
		to_file = '%s/%s/'%(target_dir,sub)

		try:
			os.makedirs(to_file)
		except:
			pass
		
		copy_cmd = '%s %s %s'%(copy_cmd_beg,from_file,to_file)
		os.system(copy_cmd)
