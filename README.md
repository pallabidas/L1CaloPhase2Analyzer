# Phase 2 L1 Calo Analyzer

## Description

   Adapted from: https://github.com/skkwan/phase2-l1Calo-analyzer.
   This repo is for running the Phase-2 calo emulator, in particular checking the digitized version of the
   calo jet emulator.

## Setup (do only once)

   ```
   cmsrel CMSSW_15_1_0
   cd CMSSW_15_1_0/src
   cmsenv
   git cms-init
   git cms-merge-topic -u pallabidas:test_GCT_TM18
   cd L1Trigger
   git clone git@github.com:pallabidas/L1CaloPhase2Analyzer.git -b 15_1_0
   scram b -j 12
   ```

## To run the emulator and create ntuples for the event display, efficiency plots

   For getting the ntuple:
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/test/
   cmsRun test-analyzer.py
   ```

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
