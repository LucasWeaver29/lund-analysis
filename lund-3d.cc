// Make jet spectra (pt, y, phi) histograms to get comfortable with fastjet and root.

#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "math.h"

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TCanvas.h"

// ROOT - interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT - saving file.
#include "TFile.h"

#include "thermalBackground.h"

using namespace Pythia8;

// Vector projection of b onto a. Vectors represented by vectors of length 3 representing x, y, z componenets
// Currently this is not being used
vector<double> vecProjection (vector<double> a, vector<double> b) {
    double abs_a = sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
    double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    
    vector<double> proj;
    for (int i = 0; i < 3; ++i) {
        proj.push_back((dot*a[i])/pow(abs_a,2));
    }
    return proj;
}

// Takes a CA clustered jet, returns a vector of kinematic 
// variables (delta, kt) vectors for each primary splitting
vector<vector<double>> decluster(fastjet::PseudoJet jet) {
    
    vector<vector<double>> kinematic_vars;

    fastjet::PseudoJet parent1;
    fastjet::PseudoJet parent2;
    // Pseudojets a and b 
    fastjet::PseudoJet* pa = &parent1;
    fastjet::PseudoJet* pb = &parent2;

    // For each jet, iteratively compare branchings
    while (jet.has_parents(parent1, parent2)) { 
        // In each case, identify the higher pt parent 
        if(parent1.pt() > parent2.pt()) {
            pa = &parent1;
            pb = &parent2;
        }
        else {
            pa = &parent2;
            pb = &parent1;
        }
        // Classify pb as the emission, and pa + pb as the emitter
        // and define kinematic variables
        
        // rap-phi distance, delta = delta_ab
        double delta = pa->delta_R(*pb);

        double kt = sin(delta) * pb->pt();

        vector<double> vars = {delta, kt};

        kinematic_vars.push_back(vars);
    
        jet = *pa;
    }
    return kinematic_vars;
}

int main(int argc, char* argv[]) {

    // number of events per bin
    // do 5000
    int nEvents = 500;
    
    TString file_name = "Lund Planes + Jet Cuts, with Background";

    // Minimum jet pT
    double pTmin = 0;
    double etaMax = 5;

    double Rparam = .4;
    
    TString xAxisName = "ln((R=" + to_string(Rparam) + "/delta)";
    TString yAxisName = "ln(kt)";

// -----------------------------------------------------------------------
    
// Preparing for background generation
    
    bool background = true;

    //std::random_device rd;
    //std::mt19937 gen(rd()); 
    std::mt19937 gen(0);    
    
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    std::normal_distribution<double> normal_dist(0.0, 1.0);

    bool flow = true;

    // For anisotropic flow distribution
    aniso_flow_distribution flow_dist;
    
    // For pT distribution
    pt_distribution pt_dist;

// -----------------------------------------------------------------------


    // dcut for exclusive jets
    // double dcut = 0;

    Pythia pythia;
    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all=on");
    pythia.readString("Next:numbershowEvent = 0");

    // Give bounds for each bin
    vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    //vector<double> ptHatMins = {30, 60, 90, 120, 150, 180};
    //vector<double> ptHatMaxs = {40, 70, 100, 130, 160, 210};
    //vector<double> ptHatMins = {0, 100};
    //vector<double> ptHatMaxs = {100, 150};
    

    // Declare 2D ROOT histograms
    TH2F *lund_akt = new TH2F("lund_akt", file_name + ", inclusive pt; " + xAxisName + "; " + yAxisName, 50, 0, 5, 50, -4, 3.5);
    
    TH3F *lund_3d = new TH3F("lund", file_name + "; " + xAxisName + "; " + yAxisName + "; jet pT", 50, 0, 5, 50, -4, 3.5, 50, 0, 200);

    // ROOT hists for pt cuts
    vector<double> ptCutMins = {10, 20, 30, 50, 70, 130, 160, 190};
    vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 160, 190, 250};
    vector<TH2F*> lundCuts;
    
    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        lundCuts.push_back(new TH2F("", file_name + ", jet pT " + ptCutMins[iCut] + " to " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName,
            50, 0, 5, 50, -4, 3.5));
    }


    // Setup fastjet analysis 
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;

    // Only interested in anti-kt initial clustering, then ca reclustering
    fastjet::JetDefinition jetDef_akt(fastjet::antikt_algorithm, Rparam, recombScheme, strat);
    
    // Make R large for reclustering CA jet definition so all jet consitutents previously identified are included in the new reclustered jet 
    fastjet::JetDefinition jetDef_ca_recluster(fastjet::cambridge_algorithm, 3, recombScheme, strat);


    // Bin loop
    for (int iBin = 0; iBin < ptHatMins.size(); ++iBin) {

        pythia.settings.readString("PhaseSpace:pTHatMin = " + to_string(ptHatMins[iBin]));
        pythia.settings.readString("PhaseSpace:pTHatMax = " + to_string(ptHatMaxs[iBin]));

        // Proton - proton collision is default
        pythia.init();

        // Event loop
        for (int n = 0; n < nEvents; ++n) {
            
            pythia.next();

            vector<fastjet::PseudoJet> fjInputs;

            // Particle loop
            for (int i=0; i<pythia.event.size(); ++i) {

                if (not pythia.event[i].isFinal()) continue;

                if (abs(pythia.event[i].eta()) > etaMax) continue;

                if (not pythia.event[i].isCharged()) continue;

                fjInputs.push_back(fastjet::PseudoJet(
                    pythia.event[i].px(), 
                    pythia.event[i].py(),
                    pythia.event[i].pz(),
                    pythia.event[i].e()));
            }


            // Weight for this bin, to account for cross section of this pT.
            double xsec = pythia.info.sigmaGen();
            //double weight = xsec / (nEvents * (ptHatMaxs[iBin]-ptHatMins[iBin]));
            double weight = xsec / nEvents;
            

            // -----------------------------------------------------------------------
            
            // Creating the background
            
            if (background) {
                // Create a vector to store all new particles
                vector<fastjet::PseudoJet> background_prtcls;

                // number of thermal particles
                // "Chosen to match multiplicities for 0-5% most central Alic Pb-Pb collisions"
                // N_ch ~ N(mu = 3922, sigma = 40)
                int N_ch = normal_dist(gen)*40 + 3922;
                //N_ch = 1000000;

                double pt, phi, eta;
                
                for (int i = 0; i < N_ch; i++) {
                    
                    pt = pt_dist.sample();

                    // Azimuthal Distribution and Elliptic Flow
                    // Guessed with n and T
                    //if (flow) phi = flow_dist.sample_swv(pt);  // Anisotropic flow according to flow-modulated distribution
                    if (flow) phi = flow_dist.sample(pt);  // Anisotropic flow according to flow-modulated distribution        
                    else phi = uniform_dist(gen)*2*M_PI;  // Uniformly 

                    // Pseudorapidity
                    // "Each particles is restricted to the ALICE central-barrel acceptance |η| < 0.9"
                    // ... So I guess pseudorapidity is assignied uniformly
                    eta = uniform_dist(gen) * 1.8 - .9;
                
                    fastjet::PseudoJet prtcl;
                    prtcl.reset_PtYPhiM(pt, eta, phi, 0);
                    
                    background_prtcls.push_back(prtcl);
                }
                
                move(background_prtcls.begin(), background_prtcls.end(), back_inserter(fjInputs));
            }
            
            // done creating background

            // -----------------------------------------------------------------------


            // Cosntruction of primary Lund plane from 
            // https://hal.science/hal-01851158/document

            fastjet::ClusterSequence clustSeq_akt(fjInputs, jetDef_akt);
            
            // Use inclusive jets here? Use pTmin?
            vector<fastjet::PseudoJet> inclusiveJets_akt = clustSeq_akt.inclusive_jets(pTmin);
            
            // Now recluster the jets identified by anti-kt according to CA, add points to lund plane
            bool cut = false;
            int iCut = 0;
            
            vector<fastjet::PseudoJet> reclusteredJets_akt;
            for (fastjet::PseudoJet jet : inclusiveJets_akt) {
                
                cut = false;

                // Find which cut this jet belongs to
                for (int i = 0; i < ptCutMins.size(); ++i) {
                    if (jet.pt() > ptCutMins[i] && jet.pt() < ptCutMaxs[i]) {
                        cut = true;
                        iCut = i;
                        break;
                    }
                }
                
                fastjet::ClusterSequence reclustered_jet(clustSeq_akt.constituents(jet), jetDef_ca_recluster);
                if (reclustered_jet.inclusive_jets().size() != 1) {
                    cout << "Warning, an anti-kt jet reclustered into " << reclustered_jet.inclusive_jets().size() << " CA jets." << endl;
                }
                
                vector<vector<double>> kin_vars = decluster(reclustered_jet.inclusive_jets()[0]);
                for (vector<double> vec : kin_vars) {
                    lund_akt->Fill(log(Rparam/vec[0]),log(vec[1]), weight);
                    lund_3d->Fill(log(Rparam/vec[0]),log(vec[1]), jet.pt(), weight);

                    if (cut) lundCuts[iCut]->Fill(log(Rparam/vec[0]),log(vec[1]), weight);


                }
            }
            
            
        }
    }

    


    TApplication app("app", nullptr, nullptr);

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
   
    c1->cd();
    lund_akt->Draw("COLZ");
    c1->Print(file_name+".pdf(","pdf");
    c1->Clear();


    for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
        c1->cd();
        lundCuts[iCut]->Draw("COLZ");
        if (iCut + 1 == ptCutMins.size()) {
            c1->Print(file_name + ".pdf)","pdf");
        } 
        else {
            c1->Print(file_name + ".pdf","pdf");
        }
        c1->Clear();
    }

    /*
    c1->cd();
    lund_3d->Draw("BOX");
    c1->Update();
    //gPad->WaitPrimitive();
    app.Run();
    */


    return 0;
}