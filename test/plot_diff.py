import argparse
import awkward as ak
import hist
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os
import uproot

plt.style.use(hep.style.CMS)

TREEPATH = "l1NtupleProducer/displayTree"

def main(input,outdir):
    if not os.path.exists(outdir):
        print(f"Creating directory {outdir}")
        os.makedirs(outdir)

    print(f"Reading in input file {input}")
    with uproot.open(input) as f:
        tree = f[TREEPATH]
    etJet = tree["gctCaloJets_et"].array()
    etDigi = tree["gctCaloJetsDigitized_etFloat"].array()
    phiJet = tree["gctCaloJets_phi"].array()
    phiDigi = tree["gctCaloJetsDigitized_phiFloat"].array()
    etaJet = tree["gctCaloJets_eta"].array()
    etaDigi = tree["gctCaloJetsDigitized_etaFloat"].array()

    #Histogram for ET diffs
    et_axis = hist.axis.Regular(12,-.03125,.03125,name="et_diff",label=r"$ET_{Float}-ET_{Digi}$")
    et_h = hist.Hist(et_axis,name="Yield")
    et_diff = etJet - etDigi
    et_h.fill(et_diff=ak.ravel(et_diff))

    #Histogram for phi diffs (note: LSB ~ 0.000767)
    phi_axis = hist.axis.Regular(8,-0.08,0.08,name="phi_diff",label=r"$\phi_{Float}-\phi_{Digi}$")
    phi_h = hist.Hist(phi_axis,name="Yield")
    phi_diff = phiJet - phiDigi
    phi_h.fill(phi_diff=ak.ravel(phi_diff*100))

    #Histogram for eta diffs (note: LSB ~ 0.000767)
    eta_axis = hist.axis.Regular(8,-0.08,0.08,name="eta_diff",label=r"$\eta_{Float}-\eta_{Digi}$")
    eta_h = hist.Hist(eta_axis,name="Yield")
    eta_diff = etaJet - etaDigi
    eta_h.fill(eta_diff=ak.ravel(eta_diff*100))

    #Plot ET hist
    fig, ax = plt.subplots(1,1,figsize=(7,4))
    hep.histplot(et_h,flow="show")
    ax.set_title(r"$E_T$ Float - Digi")
    ax.set_ylabel("Yield")
    ax.set_xlabel(r"$ET_{Float}-ET_{Digi}$")
    et_outpath = os.path.join(outdir,"ET_diff.png")
    plt.savefig(et_outpath)
    print(f"Plot saved to {et_outpath}")

    #Plot phi hist
    fig, ax = plt.subplots(1,1,figsize=(7,4))
    hep.histplot(phi_h,flow="show")
    ax.set_title(r"$\phi$ Float - Digi (x100)")
    ax.set_ylabel("Yield")
    ax.set_xlabel(r"$\phi_{Float}-\phi_{Digi}$")
    phi_outpath = os.path.join(outdir,"phi_diff.png")
    plt.savefig(phi_outpath)
    print(f"Plot saved to {phi_outpath}")

    #Plot eta hist
    fig, ax = plt.subplots(1,1,figsize=(7,4))
    hep.histplot(eta_h,flow="show")
    ax.set_title(r"$\eta$ Float - Digi (x100)")
    ax.set_ylabel("Yield")
    ax.set_xlabel(r"$\eta_{Float}-\eta_{Digi}$")
    eta_outpath = os.path.join(outdir,"eta_diff.png")
    plt.savefig(eta_outpath)
    print(f"Plot saved to {eta_outpath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot the difference between CaloJet and DigiCaloJet attributes"
    )
    parser.add_argument("--input", type=str, default="analyzer.root", help="Root file to analyze")
    parser.add_argument("--outdir", type=str, default="Outputs", help="Where to put the output files")

    args = parser.parse_args()

    main(args.input,args.outdir)