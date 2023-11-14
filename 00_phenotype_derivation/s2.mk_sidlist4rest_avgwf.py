
import os
import glob
import commands
import pandas as pd
import time
import numpy as np

start = time.time()
#seg_file = './segs/langNet31_atlas.nii.gz'
#seg_file = './segs/langNet18_CORE.nii.gz'
seg_file = './segs/AICHA.nii'

dat_dir = '/data/clusterfs/lag/users/xiakon/dat'
sid_path_list = glob.glob(os.path.join(dat_dir, ('[0-9]'*7)))
print '------------------------------------'
print len(sid_path_list), 'subject folders were found. Running...'
print '------------------------------------'
#print sid_path_list[0]

sid_list = [os.path.basename(sid) for sid in sid_path_list]

# for recording subject with avgwf, and number of regions in seg
sid_tmp = []
seg_N_tmp = []
sid_tmp2 = []

for sid in sid_list:
    print sid
    rest_clean_seg_dir = os.path.join(dat_dir, sid, 'fMRI/rfMRI.ica/seg/')
    seg_avgwf_file = os.path.join(rest_clean_seg_dir, os.path.basename(seg_file).split('.')[0]+'_rest_avgwf.txt')
    #/data/clusterfs/lag/users/xiakon/dat/5837334/fMRI/rfMRI.ica/seg/langNet31_atlas_rest_avgwf.txt
    if os.path.exists(seg_avgwf_file):
        avgwf_dat = np.loadtxt(seg_avgwf_file)
        seg_N = avgwf_dat.shape[1]
        seg_N_tmp.append(seg_N)
        sid_tmp.append(sid)
        #if seg_N ==32:
        #if seg_N ==19:
        #if seg_N ==385:
        if (seg_N ==385): #and (not os.path.exists('./mats_AICHA/fcmatzw_pos/'+sid+'_fcmatzw_pos.txt')):
            sid_tmp2.append(sid)

avgwf_log = pd.DataFrame()
avgwf_log['SID']=sid_tmp
avgwf_log['SegN'] = seg_N_tmp
avgwf_log.to_csv(os.path.basename(seg_file).split('.')[0]+'_segstats_log.csv', index=False)

sid_log = pd.DataFrame()
sid_log['SID']=sid_tmp2
sid_log.to_csv(os.path.basename(seg_file).split('.')[0]+'_segstats_sidlist.csv', index=False, header=False)

duration = (time.time() - start)
print 'Time used :', duration/60, ' mins'

