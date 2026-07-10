//#include "fastjet/UserInfo.hh"
#include "pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "TH1F.h"
#include "TCanvas.h"

#include "eventData.h"

using namespace Pythia8;

int main() {

    TString eventFileName = "events1000x9.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";

    TString output_file_name = "jet vs pc tracks";

    bool debug = false;

    if (debug) cout << "Declaring my_event_tree" << endl;
    my_event_tree met(eventFileName);
    if (debug) cout << "Declaring my_background_tree" << endl;
    my_background_tree mbt(backgroundFileName);


    double Rparam = .4;
    double part_eta_max = .9;
    double jet_eta_max = .5;
    double part_pt_min = .15;
    bool leading_jet_only = false;
    if(debug) cout << "Declaring histograms" << endl;
    // Specifically, tracks referring to background particles only
    int pt_numBins = 60;
    double pt_xMin = -.5;
    double pt_xMax = 5.6;
    TH1F* jet_tracks_pt = new TH1F("jet_tracks_pt", "Weighted pT of Background Particles in Jet Cone; Track pT; Weighted counts", pt_numBins, pt_xMin, pt_xMax);
    TH1F* pc_tracks_pt = new TH1F("pc_tracks_pt", "Weighted pT of Background Particles in PC Cone; Track pT; Weighted counts", pt_numBins, pt_xMin, pt_xMax);
    
    int area_numBins = 50;
    double area_xMin = 0;
    double area_xMax = 3;
    TH1F* akt_jet_area = new TH1F("akt_jet_area", "Area of akt jets", area_numBins, area_xMin, area_xMax);
    // assuming area of pc is piR^2
    TH1F* pc_area = new TH1F("pc_area", "Area of pc", area_numBins, area_xMin, area_xMax);


    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best);
    fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));
    fastjet::Selector jet_eta_selector = fastjet::SelectorEtaRange(-1 * jet_eta_max, jet_eta_max);
    fastjet::Selector r_selector = fastjet::SelectorCircle(Rparam);

    int numEvents = met.GetEntries();
    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }

    int counter = 1000;

    for (int iEvent = 0; iEvent < numEvents; iEvent++) {

        if (debug) cout << "iEvent " << iEvent << endl;
        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max); // iEvent, part_eta_max

        for (fastjet::PseudoJet &p : event_particles) { // index 0 for event particle, 1 for background particle
            p.set_user_index(0);
        }

        // Embedding in the background
        vector<fastjet::PseudoJet> bg_particles = mbt.get_particles(iEvent);
        for (fastjet::PseudoJet &bp : bg_particles) { // index 0 for event particle, 1 for background particle
            bp.set_user_index(1);
        }
        
        move(bg_particles.begin(), bg_particles.end(), back_inserter(event_particles));

        if (debug) cout << "Clustering event_particles" << endl;
        fastjet::ClusterSequenceArea clust_seq(event_particles, jet_def_akt, area_def);
        vector<fastjet::PseudoJet> jets = jet_eta_selector(clust_seq.inclusive_jets());
        if (leading_jet_only) jets = sorted_by_pt(jets);
        
        if (debug) cout << "Starting jet loop" << endl;
        for (fastjet::PseudoJet &jet : jets) {
            akt_jet_area->Fill(jet.area());
            
            for (fastjet::PseudoJet &particle : jet.constituents()) {
                if (particle.user_index()==1) jet_tracks_pt->Fill(particle.pt());
            }
            
            fastjet::PseudoJet pc_axis;
            pc_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));

            r_selector.set_reference(pc_axis);
            
            vector<fastjet::PseudoJet> pc_particles = r_selector(event_particles);
            for (fastjet::PseudoJet &particle : pc_particles) {
                if (particle.user_index() == 1) pc_tracks_pt->Fill(particle.pt());   
            }

            if (leading_jet_only) break;
            
        } // end jet loop

    } // end event loop
    

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    vector<TH1F*> hists_to_print = {jet_tracks_pt, pc_tracks_pt, akt_jet_area};

    for (int iHist = 0; iHist < hists_to_print.size(); iHist++) {
        c1->cd();
        hists_to_print[iHist]->Draw("HIST");
        if (iHist == 0) c1->Print(output_file_name + ".pdf(");
        else if (iHist == hists_to_print.size() - 1) c1->Print(output_file_name + ".pdf)");
        else c1->Print(output_file_name + ".pdf");
        c1->Clear();
    }

    // Now, we can go back through all the events, and look at what correction factor would account for upward background fluctuation recieved by the leading jet
    
    

}