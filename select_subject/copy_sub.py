"""
-----------------------------------------------------------------------------------------
copy_sub.py
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
cd /home/szinte/projects/retino_HCP/
python select_subject/copy_sub.py
-----------------------------------------------------------------------------------------
"""

import os
import glob
import ipdb

subs = ["192641", "105923", "111312", "926862", "182739", "167440", "789373", "690152", "148133", "467351",
        "233326", "144226", "638049", "318637", "165436", "771354", "995174", "463040", "187345", "547046",
        "150423", "167036", "104416", "536647", "397760", "552241", "169040", "973770", "116726", "130114", 
        "135124", "181636", "177746", "239136", "205220", "200614", "185442", "131722", "406836", "732243", 
        "178142", "412528", "146129", "134829", "114823", "140117", "818859", "429040", "191336", "157336", 
        "169747", "814649", "171633", "214019", "901139", "680957", "751550", "905147", "178647", "352738", 
        "195041", "971160", "191033", "330324", "927359", "320826", "654552", "158035", "115017", "249947", 
        "601127", "617748", "102311", "861456", "155938", "346137", "156334", "100610", "283543", "145834", 
        "783462", "164131", "671855", "257845", "196144", "108323", "872764", "137128", "585256", "181232", 
        "126426", "193845", "825048", "198653", "393247", "901442", "724446", "942658", "164636", "209228", 
        "214524", "263436", "177140", "389357", "102816", "572045", "859671", "118225", "782561", "833249", 
        "951457", "381038", "180533", "899885", "201515", "146432", "200210", "182436", "725751", "146937", 
        "197348", "360030", "203418", "146735", "966975", "191841", "757764", "246133", "706040", "186949", 
        "380036", "204521", "128935", "826353", "212419", "943862", "395756", "134627", "172130", "765864", 
        "541943", "159239", "898176", "175237", "385046", "192439", "176542", "581450", "878877", "878776", 
        "199655", "958976", "644246", "770352", "173334", "111514", "436845", "910241", "401422", "162935", 
        "550439", "169343", "627549", "221319", "871762", "109123", "158136", "125525", "525541", "573249", 
        "251833", "562345", "130518", "200311", "177645", "131217", "178243", "115825", "132118", "169444", 
        "365343"] 

files_to_trans = ['RETBAR_ALL_tfMRI_data_sub.nii.gz']

database_dir = '/home/raw_data/2018/visual/HCP_RETINO_HR'
target_dir = '/home/shared/2018/visual/subcortical_hcp/raw_data'

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
