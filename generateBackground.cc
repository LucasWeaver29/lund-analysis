//#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"

#include <iostream>
#include <random>

// ROOT - histogram
//#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
//#include "TLegend.h"

// ROOT - function
#include "TF1.h"

// ROOT - saving as files
#include <TFile.h>
#include <TTree.h>
//#include <TBranch.h>

#include "cmath"

using namespace std;


// Generating the thermal background laid out in 
//"Thermal Model of the Underlying Event Background in Heavy-Ion Collisions for Jet Analyses"
// and storing it as a root file to avoid regeneration

// Eyeballed from Fig 2a in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.132302
double flow_coeff(double pT = .5) {      
    if (pT < .5) {return .02;}
    if (pT < 1) {return .03;}
    if (pT < 1.5) {return .045;}
    else return .06;
}


int main() {


    //int numBackgrounds = 2700000;
    int numBackgrounds = 14000; // for a small test

    // For checking specs after generation
    TString output_name = "Background Specs - 14kb";

    // Create a RSOOT output file
    TFile file("backgrounds14kb.root", "RECREATE");
    
    // Lily restricts backgrounds to |eta| < .9. To generate a background covering a wider eta range, I use this multiplier for both eta range and particle multiplicity
    double eta_max = .9;
    bool diff_eta = eta_max == .9 ? false : true;
    double multiplier = eta_max /.9; 


    bool flow = true;

    // Create a TTree
    TTree tree("backgrounds", "Backgrounds Tree");
    tree.SetMaxTreeSize(50000000000LL);  // 50 GB per file
    
    // Define variables to store event data
    int numParticles;
    std::vector<double> pts;
    std::vector<double> phis;
    std::vector<double> etas;
    
    // Create branches for the tree
    tree.Branch("numParticles", &numParticles, "numParticles/I");
    tree.Branch("pts", &pts);
    tree.Branch("phis", &phis);
    tree.Branch("etas", &etas);

    //=================================================================
    // For checking specs after generation
    
    TH1F* nPart_hist = new TH1F("nPart_hist", "# particles per event; # of charged particles; counts", 450, 3650 * multiplier, 4100 * multiplier);
    TH1F* pt_hist = new TH1F("pt_hist", "pT Spectrum; pt; counts", 100, 0, 5);
    TH1F* phi_hist = new TH1F("phi_hist", "Azimuthal angle of UE thermal particles (anisotropic flow = true); Azimuthal angle; counts", 100, 0, 2*M_PI);
    TH1F* eta_hist = new TH1F("eta_hist", "Pseudorapidity of UE thermal particles; Pseudorapidity; counts", 100, -5, 5);


    // =================================================================
    // Prepare distrbutions

    std::random_device rd;
    std::mt19937 gen(rd()); 
    // For same 0 seed
    //std::mt19937 gen(0);    
    
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    std::normal_distribution<double> normal_dist(0.0, 1.0);


    // Anisotropic flow distribution
    // From https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.132302
    // In the paper, pT units are pT/c. But range seems to match pT genereated by pT_distribution,
    // so I'm ignoring that difference for now
    // Assumes psi = 0;
    TF1* flow_PDF = new TF1("flow_PDF", "1 + [0] * cos(2*x)", 0, 2*M_PI);
    flow_PDF->SetParName(0, "v");
    flow_PDF->SetParameter(0, .025); // If we don't feel like adjust flow coefficient for each particle pt

    
    // For pT distribution
    TF1* pt_PDF = new TF1("pt_PDF", "x * (1 + (x/([0]*[1])))^(-[0])", 0 ,5);
    pt_PDF->SetParName(0,"n");
    pt_PDF->SetParameter(0, 10);
    pt_PDF->SetParName(1,"T");
    pt_PDF->SetParameter(1, .2);
    pt_PDF->SetNpx(10000); 



    // To see the pdf
    /*
    c1->cd();
    pt_PDF->Draw();
    c1->Print("pt_PDF.pdf", "pdf");
    c1->Clear();
    */
    
    int counter = 1000;
    
    for (int iBack = 0; iBack < numBackgrounds; ++iBack) {
        
        if (iBack > counter) {
            cout << counter << " backgrounds generated" << endl;
            counter = counter + 1000;
        }

        // Clear vectors for new event
        pts.clear();
        phis.clear();
        etas.clear();

        // number of thermal particles per background
        // "Chosen to match multiplicities for 0-5% most central Alic Pb-Pb collisions"
        // N_ch ~ N(mu = 3922, sigma = 40)
        numParticles = static_cast<int>(std::round(normal_dist(gen)*40 + 3922));
        if (diff_eta) numParticles = numParticles * multiplier;
        nPart_hist->Fill(numParticles);

        double pt, phi, eta;
        // For every particle, determine pt,phi,eta
        for (int iPart = 0; iPart < numParticles; ++iPart) {
            // pt
            pt = pt_PDF->GetRandom();
            pt_hist->Fill(pt);
            pts.push_back(pt);

            // Azimuthal Distribution and Elliptic Flow
            // Guessed with n and T
            if (flow) {
                //flow_PDF->SetParameter(0, flow_coeff(pt)); // If we want to adjust distribution for flow coefficient for pt                 
                phi = flow_PDF->GetRandom();
            }  
            else phi = uniform_dist(gen)*2*M_PI;  // Uniformly 
            
            phi_hist->Fill(phi);
            phis.push_back(phi);

            // Pseudorapidity
            // "Each particles is restricted to the ALICE central-barrel acceptance |η| < 0.9"
            // ... So I guess pseudorapidity is assignied uniformly
            eta = uniform_dist(gen) * 1.8 - .9;
            if (diff_eta) eta = eta * multiplier;
            eta_hist->Fill(eta);
            etas.push_back(eta);
        }
        
       tree.Fill();

    
    }


 

    //======================================================
    //Make sure distribution looks good

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    
    c1->cd();
    nPart_hist->Draw("HIST");
    c1->Print(output_name+".pdf(","pdf");
    c1->Clear();
    

    c1->cd();
    pt_hist->Draw("HIST");
    c1->Print(output_name+".pdf(","pdf");
    c1->Clear();

    
    c1->cd();
    phi_hist->Draw("HIST");
    c1->Print(output_name+".pdf", "pdf");
    c1->Clear();
    

    c1->cd();
    eta_hist->Draw("HIST");
    c1->Print(output_name+".pdf)", "pdf");
    c1->Clear();


    // Write the tree to the file
    // Once this is done the references in the histograms can't access their data
    tree.Write();
    file.Close();


    return 0;
}