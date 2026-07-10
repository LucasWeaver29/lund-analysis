// Make jet spectra (pt, y, phi) histograms to get comfortable with fastjet and root.

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


    int nEvents = 20000;
    
    // Minimum jet pT
    double pTmin = 0;
    double etaMax = 5;

    // dcut for exclusive jets
    // double dcut = 0;

    Pythia pythia;
    
    // For using command card
    // pythia.readFile(argv[1])

    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all=on");
    //pythia.readString("PhaseSpace:pTHatMin = 40");
    pythia.readString("Next:numbershowEvent = 0");

    // Proton - proton collision is default
    pythia.init();

    // Declare pythia histograms
    //Hist pt("pT of anti-kt-clustered jets", 100, 0, 25);
    //Hist eta("Pseudorapidity of anti-kt-clustered jets", 100, -7, 7);
    // Hist phi("Azimuthal angle of anti-kt-clustered jets", 100, 0, 7);
    
    // Declare ROOT histograms
    TH1F *pt = new TH1F("pt", "Jet pT (anti-kt); pT (GeV); counts (normalized)", 100, 0, 100);
    TH1F *eta = new TH1F("eta", "Jet pseudorapidity (anti-kt); pseudorapidity ; counts (normalized)", 100, -5.5, 5.5);
    TH1F *phi = new TH1F("phi", "Jet azimuthal angle (anti-kt); azimuthal angle; counts (normalized)", 100, 0, 7);

    // Setup fastjet analysis
    fastjet::Strategy strat = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    double Rparam = .4;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme,strat);

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

        // Add each jet's stats to histograms
        for (fastjet::PseudoJet jet : inclusiveJets) {
            // ROOT histogram
            pt->Fill(jet.pt());
            eta->Fill(jet.eta());
            phi->Fill(jet.phi());
            
            // Pythia histogram
            /*
            pt.fill(jet.pt());
            eta.fill(jet.eta());
            phi.fill(jet.phi());
            */
        }
    }

    pythia.stat();

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
    c1->Print("jetSpectraOld.pdf(","pdf");
    c1->Clear();

    c1->SetLogy(0);
    c1->cd();
    eta->Draw("HIST");
    c1->Print("jetSpectraOld.pdf","pdf");
    c1->Clear();

    c1->cd();
    phi->Draw("HIST");
    c1->Print("jetSpectraOld.pdf)","pdf");
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