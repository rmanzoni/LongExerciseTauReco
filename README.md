# Tau Reconstruction exercise @ CMSPOS
In this exercise you'll see first-hand how the tau reconstruction is run on low level, Particle Flow inputs.  

To make things clearer to follow, you'll run the tau reconstruction on miniAOD samples (instead of the usual RECO/AOD). 

In this way the typical turnaround change parameters -> re-run the reconstruction -> investigate the outcome is much faster and will allow you to actually test a few different configurations and see their impact on the relevant quantities, such as efficiency and fake rejection. 

## Installation and setup
First of all you need to initialise a fresh CMSSW release and install the needed packages.

```
# setup the CMSSW release
cmsrel CMSSW_9_2_7
cd CMSSW_9_2_7/src
cmsenv

# checkout the necessary packages
git-cms-addpkg /DataFormats/PatCandidates/
git-cms-addpkg /DataFormats/TauReco/
git-cms-addpkg /PhysicsTools/PatAlgos/
git-cms-addpkg /RecoTauTag/RecoTau/

# add Jan’s repository
git remote add jan https://github.com/steggema/cmssw.git
git fetch jan

# create a new branch and make it point to Jan’s
git checkout CMSSW_9_2_X_TauRecoMiniAOD

# now compile
scram b -j16

# move to the RecoTau package
cd RecoTauTag/RecoTau/test/

# get the cfg file and run it, this will produce 'outputFULL.root' containing everything from the original files
# *plus* your new tau collections
wget https://gist.githubusercontent.com/steggema/f5b65e09ebee723f5f27b1bef53dfa03/raw/05e064b65d0fe4c3aeeb65b028ccbda857b27814/tau_miniaod.py
voms-proxy-init --voms cms
cmsRun tau_miniaod.py

# suggestion: check that the output files really contains what you want using edmDumpEventContent
```

## Produce flat root trees
Now that you have run the tau reconstruction on the original miniAOD samples and you have produced `outputFULL.root`, you need to produce handy and swift ROOT flat trees to be used afterwards for producing plots and studies.

The main piece of code is `read_taus.py`. We encourage you to go through it before running it, it is thoroughly commented and it is (should be, in case just ask!) self explanatory.

The basic idea is to run over the sample you produced, operate the geometrical matching between reco and generator-level taus and fill *two* flat ROOT trees:
* one with an entry for each generated hadronic tau
* one with an entry for each reconstructed hadronic tau

```
# clone the present package
cd $CMSSW_BASE/src/RecoTauTag/RecoTau/test/
cmsenv
git clone https://github.com/rmanzoni/TauRecoCMSPOS.git
cd TauRecoCMSPOS

# create a soft link to outputFULL.root (or alternatively move it to your convenience...)
ln -s ../outputFULL.root

# run the ntupliser! (ipython is not mandatory, but it simply happens to be more useful than the bare python in many contexts)
ipython -i read_taus.py 

# now you've finally obtained your two flat trees
tau_gen_tuple.root
tau_reco_tuple.root
```

inspect your two ntuples with ROOT's TBrowser: you'll notice that there are only a few basic quantities saved (where is the isolation?!).

Try to go through the code again, understand how to add more information to the trees and add the isolation discriminators.

# Plotting
Once you have produced the flat ntuples, you can now easily perform studies and produce plots out of them.  

For example you can inspect `plotting.py` and run it `ipython -i plotting.py`.  
It will produce two plots showing the distribution and the efficiency of reconstructed quantities with respect to the generated ones.

Use this file as an example and produce:
* pulls split by decay mode
* efficiency and fake rates for different isolation working points
* ROC curves







