#!/usr/bin/python

import json
import sys
import itertools
from multiprocessing import Pool
import subprocess

def run_job(args):
    print('starting', args)
    subprocess.call(args)
    print('done', args)

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fp:
        jobs = json.load(fp)

    params = itertools.product(jobs['sequences'], jobs['depths'],
                               jobs['inlier_thrs'], jobs['ransac_iter'],
                               (jobs['mono_left'],), (jobs['sha'],), (jobs['stereo'],))

    job_args = [('matlab', '-nodisplay', '-nosplash', '-r', "run_sequence('%s','depth_thr',%d,'inlier_thr',%d,'ransac_iter',%d, 'mono_left', %d, 'sha', '%s', 'stereo', %d); exit;" % p, '>', '/home/kreimer/KITTI/out_seq%s_depth%d_inlier%d_ransac%d_mono%d_%s_stereo%d' % p, '2>&1') for p in params]

    p = Pool(7)
    p.map(run_job, job_args)
