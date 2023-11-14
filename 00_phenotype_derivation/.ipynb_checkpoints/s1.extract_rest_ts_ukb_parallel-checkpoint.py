
# transform segs into rest space
# refs
# https://git.fmrib.ox.ac.uk/falmagro/UK_biobank_pipeline_v_1/blob/master/bb_functional_pipeline/bb_fix

import os
import glob
import commands
import time

STAGE = 'All'#'All' # '1'->'2'->'3' or 'All'

start = time.time()

#seg_file = './segs/langNet31_atlas.nii.gz'
#seg_file = './segs/langNet18_CORE.nii.gz'
seg_file = './segs/AICHA.nii'

dat_dir = '/data/clusterfs/lag/projects/lg-ukbiobank/primary_data/imaging_data/'
sid_path_list = glob.glob(os.path.join(dat_dir, ('[0-9]'*7)))
print '------------------------------------'
print len(sid_path_list), 'subject folders were found. Running...'
print '------------------------------------'

sid_list = [os.path.basename(sid) for sid in sid_path_list]

out_dir = '/data/clusterfs/lag/users/xiakon/dat'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for sid in sid_list:#[:10]:
    print sid
    rest_clean_dir = os.path.join(dat_dir, sid, 'fMRI/rfMRI.ica')
    rest_clean_file = os.path.join(rest_clean_dir, 'filtered_func_data_clean.nii.gz')
    if os.path.exists(rest_clean_file):
        rest_clean_seg_dir = os.path.join(out_dir, sid, 'fMRI/rfMRI.ica/seg/')
        if not os.path.exists(rest_clean_seg_dir):
            os.makedirs(rest_clean_seg_dir)
        
        seg_rest_file = os.path.join(rest_clean_seg_dir, os.path.basename(seg_file).split('.')[0]+'_rest.nii.gz')
        seg_avgwf_file = os.path.join(rest_clean_seg_dir, os.path.basename(seg_rest_file).split('.')[0]+'_avgwf.txt')
        warp_mni2rest_file = os.path.join(rest_clean_seg_dir, 'standard2example_func_warp')
        
        rest_ref_file = os.path.join(rest_clean_dir, 'mean_func.nii.gz')
        warp_file = os.path.join(rest_clean_dir, 'reg/example_func2standard_warp')
        if STAGE =='1':
            #cmdline = 'invwarp --ref=' + rest_ref_file + ' --warp='+warp_file + ' --out='+warp_mni2rest_file
            cmdline = 'fsl_sub -q single.q -l sgelogs invwarp --ref=' + rest_ref_file + ' --warp='+warp_file + ' --out='+warp_mni2rest_file
            #print(cmdline)
            if not os.path.exists(warp_mni2rest_file+'.nii.gz'):
                os.system(cmdline)
            # Done one is ok.
        if STAGE =='2':
            #cmdline1 = 'applywarp --ref=' + rest_ref_file + ' --in=' + seg_file + ' --out=' + seg_rest_file + ' --warp='+warp_mni2rest_file+' --interp=nn'
            #cmdline1 = 'fsl_sub -q single.q -l sgelogs applywarp --ref=' + rest_ref_file + ' --in=' + seg_file + ' --out=' + seg_rest_file + ' --warp='+warp_mni2rest_file+' --interp=nn'
            cmdline1 = 'qsub -q single.q -l sgelogs applywarp -F \"--ref=' + rest_ref_file + ' --in=' + seg_file + ' --out=' + seg_rest_file + ' --warp='+warp_mni2rest_file+' --interp=nn\"'
            #print(cmdlinei1)
            if not os.path.exists(seg_rest_file):
                os.system(cmdline1)
        if STAGE =='3':
            #cmdline2 = 'mri_segstats --i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file
            #cmdline2 = 'fsl_sub -q single.q -l sgelogs mri_segstats --i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file
            cmdline2 = 'qsub -q single.q -l sgelogs mri_segstats -F \"--i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file + '\"'
            #print(cmdline2)
            if not os.path.exists(seg_avgwf_file):
                os.system(cmdline2)
        if STAGE=='All':
            cmdline = 'invwarp --ref=' + rest_ref_file + ' --warp='+warp_file + ' --out='+warp_mni2rest_file
            if not os.path.exists(warp_mni2rest_file+'.nii.gz'):
                os.system(cmdline) # Done once is ok
            cmdline1 = 'applywarp --ref=' + rest_ref_file + ' --in=' + seg_file + ' --out=' + seg_rest_file + ' --warp='+warp_mni2rest_file+' --interp=nn'
            if not os.path.exists(seg_rest_file):
                os.system(cmdline1)
            cmdline2 = 'fsl_sub -q single.q -l sgelogs mri_segstats --i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file
            #cmdline2 = 'qsub -q single.q -l sgelogs mri_segstats -F \"--i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file + '\"'
            #cmdline2 = 'mri_segstats --i '+ rest_clean_file + ' --seg ' + seg_rest_file +'  --avgwf ' +seg_avgwf_file
            if not os.path.exists(seg_avgwf_file):
                os.system(cmdline2)

duration = (time.time() - start)
print 'Time used :', duration/60, ' mins'

