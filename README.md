# Phase 2 L1 Calo Analyzer

## Description

   Adapted from: https://github.com/skkwan/phase2-l1Calo-analyzer
   This repo is for running the Phase-2 calo emulator.

## Setup (do only once)

   ```
   cmsrel CMSSW_12_3_0_pre4
   cd CMSSW_12_3_0_pre4/src
   cmsenv
   git cms-init
   git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v3.4.44
   git cms-merge-topic pallabidas:test-pfclusters-CMSSW_12_3_0_pre4
   scram b -j 12
   ```

   Inside `L1Trigger` directory:
   ```
   git clone git@github.com:pallabidas/L1CaloPhase2Analyzer.git
   scram b -j 12
   ```

## To run the emulator and create ntuples for the event display, efficiency plots

   For running the emulator:
   ```
   cd L1Trigger/L1CaloTrigger/test/
   cmsRun test_Phase2L1CaloEGammaEmulator.py
   ```

   To call the firmware-based emulator and make an n-tuple for event display:
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/test/
   cmsRun test-l1tEventDisplayGenerator.py   
   ```

   Use the `L1EventDisplay.root` file as input inside `figures/eventDisplays/plotEventDisplayPhaseIIPFclusters.C`
   Get plot by running `./makedisplay.sh` after suitably changing save location inside the script.


   For the efficiency plots:
   ```
   cd L1Trigger/L1CaloPhase2Analyzer/test/
   cmsRun test-analyzer.py
   ```

   This produces `analyzer.root`.
   Use this file as input in `figures/efficiencyPlots/makeEfficienciesPlotPF.cpp`.
   Get plots using `root -l -b -q makeEfficienciesPlotPF.cpp` after suitably changing save location inside the script.

