#!/usr/bin/env python
import random
#tot_job1=13
#tot_job2=8
tot_job1=7
tot_job2=61
for jobid1 in range(tot_job1):
    for jobid2 in range(tot_job2):
        jobname = "scan_profile_llh_{:02d}_{:02d}".format(jobid1,jobid2)
        print("JOB {} /data/user/hhamdaoui/fit_pass2_sys/submission_profile2d/OneJob.submit".format(jobname))
        print('VARS {} JOBID1="{}" TJOB1="{}" JOBID2="{}" TJOB2="{}" OUT="astronormindex"'.format(jobname, jobid1, tot_job1, jobid2, tot_job2))
