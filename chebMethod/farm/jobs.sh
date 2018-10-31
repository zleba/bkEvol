#condor_submit -name DGLAP     tag=DGLAP    nJobs=100  jobs.submit
#condor_submit -name DGLAPadd  tag=DGLAPadd nJobs=150  jobs.submit
#condor_submit -name fullBFKL  tag=fullBFKL nJobs=600  jobs.submit

condor_submit -name simpleBFKL  tag=simpleBFKL  nJobs=100  jobs.submit
