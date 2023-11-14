
import os
import glob
import commands
import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt

start = time.time()
#seg_file = './segs/langNet31_atlas.nii.gz'
#seg_file = './segs/langNet18_CORE.nii.gz'
seg_file = './segs/AICHA.nii'
#mat_dir = '/home/xiakon/MPI_workspace/lg-ukbiobank/analysis/xiangzhen/LangNet/mats'
mat_dir = '/home/xiakon/MPI_workspace/lg-ukbiobank/analysis/xiangzhen/FuncNet_AICHA/mats_'+os.path.basename(seg_file).split('.')[0]
subdir = 'fcmatzw_pos'

#sid_list_file = 'langNet31_atlas_segstats_sidlist.csv'
sid_list_file = os.path.basename(seg_file).split('.')[0]+'_segstats_sidlist.csv'
sid_list = [sid.strip() for sid in open(sid_list_file)]

N_sid = len(sid_list)
print '------------------------------------'
print N_sid, 'subjects\' mats were found. Running...'
print '------------------------------------'


out_dir = os.path.join(mat_dir, 'avgs')
print out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
out_file_str = 'avgs_'+subdir

N_counter = 0
for i, sid in enumerate(sid_list):
    #print sid
    mat_file = os.path.join(mat_dir, subdir, sid+'_'+subdir+'.txt')
    if os.path.exists(mat_file):
        mat_dat = np.loadtxt(mat_file)
        if not np.isnan(mat_dat).any():
            N_counter+=1
            np.fill_diagonal(mat_dat, np.nan)
            if i==0:
                mat_tmp = mat_dat
            else:
                mat_tmp+=mat_dat
        else:
            print sid

avg_mat = mat_tmp/N_counter

# save
np.savetxt(os.path.join(out_dir, out_file_str+'_N'+str(N_counter)+'.txt'), avg_mat)

# save plot
#plt.imshow(avg_mat, interpolation='nearest', origin='lower', cmap='bwr',vmin=0, vmax=1)
plt.imshow(avg_mat*(avg_mat>0), interpolation='nearest', origin='lower', cmap='bwr',vmin=0, vmax=1)
plt.colorbar()
##plt.show()
plt.savefig(os.path.join(out_dir, out_file_str+'_N'+str(N_counter)+'.png'), dpi=300, bbox_inches='tight')


duration = (time.time() - start)
print 'Time used :', duration/60, ' mins'


