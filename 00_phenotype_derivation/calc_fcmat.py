#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Docstrings for the modula.

"""
import argparse
import numpy as np
import os
#from nkpi.base import readsessid

def readsessid(fsessid):
    """Read session indetifier file.
    input:
    fsessid:  file for session indentifier
    output:
    sessid: a list for session indentifier 
    """
    fsessid = open(fsessid)
    sessid  = [line.strip() for line in fsessid]
    return sessid

def main():
    parser = argparse.ArgumentParser(description = 'calc_fcmat.',
                                     prog = 'calc_fcmat')
    parser.add_argument('-sd',
                        dest = 'sessdir',
                        required = True,
                        metavar = 'sessid-dir',
                        help = 'session data dir.')
    parser.add_argument('-sf',
                        dest = 'sessfile',
                        required = True,
                        metavar = 'sessid-file',
                        help = 'an input file containing subject id list.')
    parser.add_argument('-in',
                        dest = 'infile',
                        required = True,
                        metavar = 'infile',
                        help = 'infile.')
    parser.add_argument('-outdir',
                        dest = 'outdir',
                        required = True,
                        metavar = 'outdir',
                        help = 'out-dir.')

    args = parser.parse_args()
    
    
    sesslist = readsessid(args.sessfile)
    subdirs = ['fcmat', 'fcmatz', 'fcmatw', 'fcmatw_pos', 'fcmatzw_pos']
    for subdir in subdirs:
        tmpdir = os.path.join(args.outdir, subdir)
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    
    #
    #roi_list = [0,1,3,4,5,6,7,8,11,13,14,15,16,23,26,32,33,34,35,36,37,38,39]
    #

    for sess in sesslist:
        print sess
        tsfile = os.path.join(args.sessdir, sess, args.infile)
        tsdat = np.loadtxt(tsfile)
        #
        #tsdat = tsdat[:,roi_list]
        tsdat = tsdat[:,1:] # remove meants for segid=0
        #
        fcmat = np.corrcoef(tsdat.T)
        np.fill_diagonal(fcmat, 0.99999) # Set diagonal to 0.99999 to avoid divide 0 issue, and set to np.nan later
        outfile = os.path.join(args.outdir, subdirs[0], sess + '_fcmat.txt')
        np.savetxt(outfile, fcmat, fmt='%.4e')

        fcmatz = 0.5*np.log((1.0+fcmat)/(1.0-fcmat))
        #fcmatz[fcmatz==np.inf] = fcmatz[fcmatz<np.inf].max()
        outfile = os.path.join(args.outdir, subdirs[1], sess + '_fcmatz.txt')
        np.savetxt(outfile, fcmatz, fmt='%.4e')
        
        fcmatw = (fcmat + 1.0)/2.0 # For weighted network
        outfile = os.path.join(args.outdir, subdirs[2], sess + '_fcmatw.txt')
        np.savetxt(outfile, fcmatw, fmt='%.4e')
        
        fcmat[fcmat<0]=0 # For weighted network by remove negative correlations
        outfile = os.path.join(args.outdir, subdirs[3], sess + '_fcmatw_pos.txt')
        np.savetxt(outfile, fcmat, fmt='%.4e')
        
        fcmatz[fcmatz<0] = 0 # For weighted network by remove negative correlations
        outfile = os.path.join(args.outdir, subdirs[4], sess + '_fcmatzw_pos.txt')
        np.savetxt(outfile, fcmatz, fmt='%.4e')

if __name__ == '__main__':
    main()
