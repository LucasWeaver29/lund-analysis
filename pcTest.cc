#include <fastjet/UserInfo.hh>

#include "backSubTools.h"

int main() {

    TString eventFileName = "events1000x9, bin weight.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";

    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);


    double Rparam = .4;
    double part_eta_max = .9;
    double jet_eta_max = .5;
    double jet_pt_min = 20;
    double jet_pt_max = 200;
    double part_pt_min = .15;

    bool leading_jet_only = false;

    // Specifically, tracks referring to background particles only
    int pt_numBins = 50;
    double pt_xMin = 0;
    double pt_xMax = 50;
    // These look at, in general, how the pt of constituents of jet cone and pc compare. Really, we have to compare on a jet pt basis, or on a bin basis. 
    TH1F* jet_tracks_pt = new TH1F("jet_tracks_pt", "Weighted pT of Background Particles in Jet Cone; Track pT; Weighted counts", pt_numBins, pt_xMin, pt_xMax);
    TH1F* pc_tracks_pt = new TH1F("pc_tracks_pt", "Weighted pT of Background Particles in PC Cone; Track pT; Weighted counts", pt_numBins, pt_xMin, pt_xMax);
    
    int area_numBins = 50;
    double area_xMin = 0;
    double area_xMax = 20;
    TH1F* akt_jet_area = new TH1F("akt_jet_area", "Area of akt jets", area_numBins, area_xMin, area_xMax);
    //TH1F* pc_area = new TH1F("pc_area", "Area of pc", area_numBins, area_xMin, area_xMax);


    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best);
    fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));
    fastjet::Selector jet_eta_selector = fastjet::SelectorEtaRange(-1 * jet_eta_max, jet_eta_max);
    
    my_pc_subtractor pc_subtractor(Rparam);

    int numEvents = met.GetEntries();
    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }

    int counter = 1000;

    for (int iEvent = 0; iEvent < numEvents; iEvent++) {

        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }


        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max, part_pt_min); // iEvent, part_eta_max

        for (fastjet::PseudoJet &p : event_particles) { // index 0 for event particle, 1 for background particle
            p.set_user_index(0);
        }

        // Embedding in the background
        vector<fastjet::PseudoJet> bg_particles = mbt.get_particles(iEvent, part_pt_min);
        for (fastjet::PseudoJet &bp : bg_particles) { // index 0 for event particle, 1 for background particle
            bp.set_user_index(1);
        }

        
        move(bg_particles.begin(), bg_particles.end(), back_inserter(event_particles));


        fastjet::ClusterSequenceArea clust_seq(event_particles, jet_def_akt, area_def);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(jet_eta_selector(clust_seq.inclusive_jets()));
        
        for (fastjet::PseudoJet& jet : jets) {

            akt_jet_area->Fill(jet.area());

            for (fastjet::PseudoJet& particle : jet.constituents) {
                if (particle.get_user_index() == 1) jet_tracks_pt->Fill(particle.pt(), met.bin_weight);
            }

            vector<fastjet::PseduoJet> pc_particles = pc_subtractor.get_pc_particles(jet, event_particles);
            for (fastjet::PseudoJet& particle : pc_particles) {
                if (particle.get_user_index() ==  1) pc_tracks_pt->Fill(particle.pt(), met.bin_weight);
            }

            if (leading_jet_only) break;
        }



    }

}