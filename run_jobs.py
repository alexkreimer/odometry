#!/usr/bin/python
from __future__ import print_function
import json
import sys
import itertools
from multiprocessing import Pool
import subprocess
import uuid

def run_job(args):
    print('starting', ' '.join(args))
    subprocess.call(args)
    print('done', ' '.join(args))

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fp:
        data = json.load(fp)

    cmd = "run_sequence(%s); exit;"
    args = []
    for jobs in data:
        opts = jobs['options']
        opt_list = []
        for opt, val in opts.iteritems():
            if type(val) == int:
                opt_list.append("'%s',%d" % (opt, val))
            else:
                opt_list.append("'%s','%s'" % (opt, str(val)))
        
        opt_str = ','.join(opt_list)

        opt_prod = itertools.product(jobs['sequences'], jobs['depths'], jobs['inlier_thrs'],
                                     jobs['ransac_iter'])
        opt_list = ["'%s','depth_thr',%d,'inlier_thr',%d,'ransac_iter',%d" % opt for opt in opt_prod]

        final_opts = ['%s,%s' % (opt, opt_str) for opt in opt_list]

        for opt in final_opts:
            out_file = uuid.uuid4()
            args.append(('matlab', '-nodisplay', '-nosplash', '-r', cmd % opt , '>',
                     '/home/kreimer/KITTI/out_%s' % out_file, '2>&1'))

    p = Pool(7)
    p.map(run_job, args)
