// ROOT - files
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <iostream>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"


int main() {

    TString root_file_name = "jet energy plots 9k, leading jet only.root";
    //TString root_file_name = "Lund, bg, area sub - 45000.root";

    TString output_name = "jet energy plots 9k, leading jet only, resized";
    //TString output_name = "Lund, bg, area sub - 45000, resized";

    // Open the event ROOT file
    TFile rootFile(root_file_name, "READ");
    
    if (rootFile.IsZombie()) {
        std::cout << "Error: Could not open ROOT file" << std::endl;
        return 1;
    }

    // tree name for jes_jer: jet_energy_tree
    // tree name for makeLund: lund_tree
    // Get the TTree, INSERT TTREE NAME HERE
    TTree* tree = (TTree*)rootFile.Get("jet_energy_tree");
    if (!tree) {
        std::cout << "Error: Could not find tree in ROOT file" << std::endl;
        return 1;
    }

    // For jes_jer
    TH2F* hist = nullptr;
    TGraph* jes = nullptr;
    TGraph* jer = nullptr;
    

    // for just lund
    //TH2F* hist = nullptr;

    // jes_jer branch name: "2dhist"
    // makeLund branch name: "lund"
    tree->SetBranchAddress("2dhist", &hist);
    //tree->SetBranchAddress("jes", &jes);
    //tree->SetBranchAddress("jer", &jer);

    tree->GetEntry(0);

    //hist->Rebin2D(2); // combine 2 bins into 1

    // Make new JES, JER spectra from 2dhist
    jes = new TGraph();
    jes->SetTitle("Jet Energy Scale vs pt True; pt True; Jet energy scale (GeV)");
    jer = new TGraph();
    jer->SetTitle("Jet Energy Resolution vs pt True; pt True; Jet energy resolution (GeV)");
    
    TH1D* xProj = hist->ProjectionX("xProj");
    std::cout << "X Underflow bin contents: " << xProj->GetBinContent(0) << std::endl;
    std::cout << "X Overflow bin contents: " << xProj->GetBinContent(xProj->GetNbinsX() + 1) << std::endl;

    TH1D* yProj = hist->ProjectionY("yProj");
    std::cout << "Y Underflow bin contents: " << yProj->GetBinContent(0) << std::endl;
    std::cout << "Y Overflow bin contents: " << yProj->GetBinContent(yProj->GetNbinsX() + 1) << std::endl;

    //Zoom in
    hist->GetXaxis()->SetRangeUser(10, 150);
    hist->GetYaxis()->SetRangeUser(-1,3);

    double xMax = 200;
    double xMin = 0;
    double nBinsX = 400;
    double binWidth = (xMax - xMin) / nBinsX;

    // Since I zoomed in
    xMin = 10;
    
    // Loop through 2D hist to extract scale and resolution
    for (int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin) {
        // Use weight for jes, jer vs pt true? -> Don't think so, because weight was accounted for in creation of 2D sub_vs_true
        TH1D* slice = hist->ProjectionY("slice", iBin, iBin);
        if (slice->GetEntries() == 0) continue;
        // plots point in the middle of the bin
        jes->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetMean());
        jer->AddPoint(iBin*binWidth + xMin + binWidth*.5, slice->GetStdDev());
    }
    


    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    
    c1->cd();
    //c1->SetLogz(1);
    hist->Draw("COLZ");
    c1->Print(output_name+".pdf(","pdf");
    c1->Clear();

    
    c1->cd();
    jes->Draw();
    c1->Print(output_name+".pdf","pdf");
    c1->Clear();

    c1->cd();
    jer->Draw();
    c1->Print(output_name+".pdf)","pdf");
    c1->Clear();

    


    return 0;
}
