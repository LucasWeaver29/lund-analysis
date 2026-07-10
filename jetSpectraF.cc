// Make jet spectra (pt, y, phi) histograms to get comfortable with fastjet and root.
// Use ptHatMin, Max to get data for many different pts,
// then combine based on each bins relative cross section

// Does all of this from event data stored in Root file

#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// ROOT - files
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <iostream>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

// ROOT - interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT - saving file.
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

    TString output_name = "jetSpectra - 1000x9pTHat.pdf";
    
    // Get bin specifications here from storeWithRoot
    // Make sure recreation of iBin later down also matches
    vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    //vector<double> ptHatMins = {10, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {50, 70, 90, 110, 130, 150, 250};

    // Minimum jet pT
    double pTmin = 0;
    //double etaMax = 5;


    //==================================================================
    // Preparing to read from ROOT file

    // Open the ROOT file
    TFile rootFile("events1000x9_ptHat.root", "READ");
    if (rootFile.IsZombie()) {
        std::cout << "Error: Could not open ROOT file" << std::endl;
        return 1;
    }

    // Get the TTree
    TTree* tree = (TTree*)rootFile.Get("events");
    if (!tree) {
        std::cout << "Error: Could not find 'events' tree in ROOT file" << std::endl;
        return 1;
    }

    // Define variables to store event data
    double ptHatMin, ptHatMax, crossSection, weight;
    double pTHat;
    int numParticles;
    std::vector<double>* px = nullptr;
    std::vector<double>* py = nullptr;
    std::vector<double>* pz = nullptr;
    std::vector<double>* energy = nullptr;
    std::vector<int>* pdgId = nullptr;

    tree->SetBranchAddress("ptHatMin", &ptHatMin);
    tree->SetBranchAddress("ptHatMax", &ptHatMax);
    tree->SetBranchAddress("crossSection", &crossSection);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("numParticles", &numParticles);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("pdgId", &pdgId);
    tree->SetBranchAddress("pTHat", &pTHat);

    // ======================================================
    // Prepare ROOT analysis
    
    // Declare ROOT histograms
    TH1F *pt = new TH1F("pt", "Jet pT (anti-kt) - Weights; pT (GeV); counts", 400, 0, 200);
    TH1F *pt_noWeights = new TH1F("pt", "Jet pT (anti-kt) - No weights; pT (GeV); counts", 400, 0, 200);
    TH1F *eta = new TH1F("eta", "Jet pseudorapidity (anti-kt); pseudorapidity ; counts", 100, -5.5, 5.5);
    TH1F *phi = new TH1F("phi", "Jet azimuthal angle (anti-kt); azimuthal angle; counts", 100, 0, 7);
    TH1F *weights = new TH1F("weights", "Event Weights; event weights; counts", 300, 0, .0025);
    //TH1F *weights_zoomed = new TH1F("weights_zoomed", "Event Weights, Zoomed in; event weights; counts", 300, 0, .000001);

    TH1F *jetpt_vs_pTHat = new TH1F("jetpt_vs_pTHat", "jetpt_vs_pTHat", 100, -5, 5);


    // ROOT histogram for each bin
    vector<TH1F*> binHists = {};
    for (int i = 0; i < ptHatMins.size(); ++i) {
        TString num = "" + i;
        TString name = "pT bin range: " + to_string(ptHatMins[i]) + "to" + to_string(ptHatMaxs[i]);
        binHists.push_back(new TH1F(name, name,400, 0, 200));
    }
    
    vector<TH1F*> binHists_noW = {};
    for (int i = 0; i < ptHatMins.size(); ++i) {
        TString num = "" + i;
        TString name = "pT bin range (unweighted): " + to_string(ptHatMins[i]) + "to" + to_string(ptHatMaxs[i]);
        binHists_noW.push_back(new TH1F(name, name,400, 0, 200));
    }
    // ======================================================


    // Setup fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double Rparam = .4;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme,strat);

    int numEvents = tree->GetEntries();
    
    int counter = 1000;
    
    int iBin = 0;
    double currentBin = ptHatMins[0];


    // Event loop for events stored in ROOT file
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        tree->GetEntry(iEvent);
        weights->Fill(weight);
        //weights_zoomed->Fill(weight);
        if (weight == 0) cout << "An event weight was 0!" << endl;

        if (ptHatMin>currentBin) {
            ++iBin;
            currentBin = ptHatMin;
        }

        vector<fastjet::PseudoJet> fjInputs;

        // Particle loop. These are the final, charged, |etc| < max particles saved by storeWithRoot
        for (int iPart=0; iPart<numParticles; ++iPart) {

            fjInputs.push_back(fastjet::PseudoJet((*px)[iPart], (*py)[iPart],(*pz)[iPart],(*energy)[iPart]));
        }

        // Fastjet analysis
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        
        vector<fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets(pTmin);
        // vector<fastket::PseudoJet> exclusiveJets = clustSeq.exclusive_jets(dcut);
        

        // Add each jet's stats to histograms
        for (fastjet::PseudoJet jet : inclusiveJets) {
            // Try to take out that one outlier
            /*
            if (iBin == 0 && jet.pt() > 100) {
                continue;
            }
            */
            jetpt_vs_pTHat->Fill(jet.pt()/pTHat);
            
            // ROOT histogram
            pt->Fill(jet.pt(), weight);
            pt_noWeights->Fill(jet.pt());
            eta->Fill(jet.eta(), weight);
            phi->Fill(jet.phi(), weight);
            
            binHists[iBin]->Fill(jet.pt(), weight);
            binHists_noW[iBin]->Fill(jet.pt());
        }
    }
    

    // pythia.stat();

    // Normalizing histograms
    /*
    pt->Scale(1.0/pt->Integral());
    eta->Scale(1.0/eta->Integral());
    phi->Scale(1.0/phi->Integral());
    */

    //pt->SetLogy(1);
    //eta->SetLogy(1);
    //phi->SetLogy(1);

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    
    
    // Weighted pt
    c1->SetLogy(1);
    c1->cd();
    pt->Draw("HIST" "PLC");
    
    // Looping through histograms for each weighted bin
    TLegend* legend = new TLegend(.1, .7, .28, .9);
    legend->AddEntry(pt, "Total", "l");
    for (int i = 0; i < binHists.size(); ++i) {
        binHists[i]->Draw("HIST" "SAME" "PLC");
        TString caption = "" + to_string((int)ptHatMins[i]) + " to " + to_string((int)ptHatMaxs[i]);
        legend->AddEntry(binHists[i], caption, "l");
    }
    legend->Draw();
    c1->Print(output_name + "(","pdf");
    c1->Clear();


    // Unweighted pt
    c1->SetLogy(1);
    c1->cd();
    pt_noWeights->Draw("HIST" "PLC");

    // Looping through histograms for each unweighted bin
    TLegend* legend2 = new TLegend(.1, .7, .28, .9);
    legend2->AddEntry(pt, "Total", "l");
    for (int i = 0; i < binHists_noW.size(); ++i) {
        binHists_noW[i]->Draw("HIST" "SAME" "PLC");
        TString caption = "" + to_string((int)ptHatMins[i]) + " to " + to_string((int)ptHatMaxs[i]);
        legend2->AddEntry(binHists[i], caption, "l");
    }
    legend2->Draw();
    c1->Print(output_name, "pdf");
    c1->Clear();

    c1->SetLogy(0);
    c1->cd();
    eta->Draw("HIST");
    c1->Print(output_name,"pdf");
    c1->Clear();

    c1->cd();
    phi->Draw("HIST");
    c1->Print(output_name,"pdf");
    c1->Clear();

    /*
    c1->cd();
    weights_zoomed->Draw("HIST");
    c1->Print(output_name,"pdf");
    c1->Clear();
    */

    c1->cd();
    jetpt_vs_pTHat->Draw("HIST");
    c1->Print(output_name,"pdf");
    c1->Clear();

    c1->cd();
    weights->Draw("HIST");
    c1->Print(output_name+")","pdf");
    c1->Clear();

    
    rootFile.Close();
    
    return 0;
}