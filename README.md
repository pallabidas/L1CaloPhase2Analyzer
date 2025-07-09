# Phase 2 L1 Calo Analyzer

## Description

   Forked from https://github.com/pallabidas/L1CaloPhase2Analyzer.git, branch 13_3_0_calojet, which in turn
   is adapted from: https://github.com/skkwan/phase2-l1Calo-analyzer.
   This repo is for running the Phase-2 calo emulator, in particular checking the digitized version of the
   calo jet emulator.

## Prerequisites to Running

   CMSSW_15_0_0_pre3 is meant to be run with an el8 OS. If you are running on a non-el8 machine, you can
   use an el8 container. For example, one can run the following command from the project root directory (the
   same directory as the one in which you will run `cmsrel`):
   ```
   export IMAGE=cms:rhel8-x86_64

   singularity exec -B ${PWD}:/srv -B /etc/condor -B /etc/grid-security/certificates \
   -B /cvmfs -B /hdfs -B /scratch -B /afs --pwd /srv /cvmfs/singularity.opensciencegrid.org/cmssw/${IMAGE} \
   /bin/bash --rcfile /srv/.bashrc
   ```
   where your current directory also contains a file called `.bashrc` with the following contents:
   ```
   unset SCRAM_ARCH
   source /cvmfs/cms.cern.ch/cmsset_default.sh
   ```
   Then, before proceeding, run `./shell` in the command line. This must be done every time before using the
   code in this repository. Note that you can run `exit` to exit the container.

## Setup (do only once)

   ```
   cmsrel CMSSW_15_0_0_pre3
   cd CMSSW_15_0_0_pre3/src
   cmsenv
   git cms-init
   git cms-addpkg L1Trigger/L1TCalorimeter
   git cms-addpkg DataFormats/L1TCalorimeterPhase2
   cd L1Trigger
   git clone https://github.com/rpsimeon34/L1CaloPhase2Analyzer.git -b 15_0_0_pre3_calojet
   cd ../
   git remote add digi-repo https://github.com/rpsimeon34/cmssw.git
   git fetch digi-repo
   git checkout -b digi-branch --track digi-repo/from-CMSSW_15_0_0_pre3-QuickDigi
   git checkout digi-branch
   git pull
   scram b -j 12
   ```

## To run the emulator and create ntuples for the event display, efficiency plots

   For getting the ntuple:
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/test/
   cmsRun test-analyzer.py
   ```
   To generate plots showing the difference between digitized and non-digitized jet attributes, run the
   following. Note that this script requires packages not included by default in the above-named container
   image:
   ```
   python3 plot_diff.py
   ```

   The remainder of this README.md is leftover from the source repository - it is not guaranteed to work here.

   For the GCT jet efficiency plots using the ntuple as input (changing file paths needed in plotting script):
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/figures/efficiencyPlots/
   root -l -b -q makeEfficienciesPlotJet.cpp
   ```

   For the event display plots using the ntuple as input (changing file paths needed in plotting script):
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/figures/eventDisplays/
   root -l -b -q plotEventDisplayPhaseIICaloJets.C
   ```
