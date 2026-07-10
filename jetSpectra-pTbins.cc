// Make jet spectra (pt, y, phi) histograms to get comfortable with fastjet and root.
// Use ptHatMin, Max to get data for many different pts,
// then combine based on each bins relative cross section

#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"

// ROOT - interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT - saving file.
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

    // Events per bin
    int nEvents = 4000;
    
    // Minimum jet pT
    double pTmin = 0;
    double etaMax = 5;


    Pythia pythia;
    
    // int binSize = 10;
    
    // Give bounds for each bin
    vector<double> ptHatMins = {30, 60, 90, 120, 150};
    vector<double> ptHatMaxs = {40, 70, 100, 130, 160};

    // Declare ROOT histograms
    TH1F *pt = new TH1F("pt", "Jet pT (anti-kt); pT (GeV); counts (normalized)", 100, 0, 200);
    TH1F *eta = new TH1F("eta", "Jet pseudorapidity (anti-kt); pseudorapidity ; counts (normalized)", 100, -5.5, 5.5);
    TH1F *phi = new TH1F("phi", "Jet azimuthal angle (anti-kt); azimuthal angle; counts (normalized)", 100, 0, 7);

    // Setup fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double Rparam = .4;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme,strat);


    // For using command card
    // pythia.readFile(argv[1])

    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all=on");
    pythia.readString("Next:numbershowEvent = 0");
    //pythia.readString("PhaseSpace:pTHatMin = 70");
    //pythia.readString("PhaseSpace:pTHatMax = 75");
 
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

                if(not pythia.event[i].isCharged()) continue;

                fjInputs.push_back(fastjet::PseudoJet(pythia.event[i].px(), 
                    pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
            }

            // Fastjet analysis
            fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
            
            vector<fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets(pTmin);
            // vector<fastket::PseudoJet> exclusiveJets = clustSeq.exclusive_jets(dcut);

            // Weight for this bin, to account for cross section of this pT.
            double xsec = pythia.info.sigmaGen();
            double weight = xsec / (nEvents * (ptHatMaxs[iBin]-ptHatMins[iBin]));

            // Add each jet's stats to histograms
            for (fastjet::PseudoJet jet : inclusiveJets) {
                // ROOT histogram
                pt->Fill(jet.pt(), weight);
                eta->Fill(jet.eta(), weight);
                phi->Fill(jet.phi(), weight);
                
                // Pythia histogram
                /*
                pt.fill(jet.pt());
                eta.fill(jet.eta());
                phi.fill(jet.phi());
                */
            }
        }
    }

    // pythia.stat();

   // Normalizing histograms
    pt->Scale(1.0/pt->Integral());
    eta->Scale(1.0/eta->Integral());
    phi->Scale(1.0/phi->Integral());

    //pt->SetLogy(1);
    //eta->SetLogy(1);
    //phi->SetLogy(1);


    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetLogy(1);
    
    // Saving each histogram manually
    c1->cd();
    pt->Draw("HIST");
    c1->Print("jetSpectra-pTbins.pdf(","pdf");
    c1->Clear();

    c1->SetLogy(0);
    c1->cd();
    eta->Draw("HIST");
    c1->Print("jetSpectra-pTbins.pdf","pdf");
    c1->Clear();

    c1->cd();
    phi->Draw("HIST");
    c1->Print("jetSpectra-pTbins.pdf)","pdf");
    c1->Clear();
    

    // A loop to save many histograms
    /*
    vector<TH1F*> hists = {pt, eta, phi};

    char file_name[] = "jetSpectra.pdf";
    
    char file_open[strlen(file_name)+2];
    strcpy(file_open, file_name);
    strcat(file_open,"(");

    char file_close[strlen(file_name)+2];
    strcpy(file_close, file_name);
    strcat(file_close, ")");

    for(int i = 0; i < hists.size(); ++i ) {        
        c1->cd();

        hists[i]->Draw("HIST");
        
        if (i==0) c1->Print(file_open,"pdf");
        else if ((i+1) == hists.size()) c1->Print(file_close,"pdf");
        else c1->Print(file_name, "pdf");
        
        c1->Clear();
    }
    
    // clean up
    for(int i = 0 ; i<hists.size(); ++i) {
        delete hists[i];
    }
    */
    /*
    // Showing ROOT histograms
    pt->Draw();
    eta->Draw();
    phi->Draw();

    gPad->WaitPrimitive();
    */
    
    // cout << pt << eta << phi;

    return 0;
}