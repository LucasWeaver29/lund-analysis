

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"


#include "TH1F.h"
#include "TCanvas.h"

#include "eventData.h"

using namespace std;

int main() {


    TString event_file_name = "event_test.root";
    vector<TString> background_file_names = {"backgrounds14ka.root", "backgrounds14kb.root"};

    TString output_name = "Compare Area Methods";

   
    bool debug = false;

    double Rparam = .4;
    double part_eta_max = .9;
    double jet_eta_max = .5;


    int area_nBins = 100;
    int area_xMin = 0;
    int area_xMax = 4;
    TH1F* areas_by_event = new TH1F("area_by_event", "Jet Areas When Found Event-wide", area_nBins, area_xMin, area_xMax);
    TH1F* areas_by_akt_jet = new TH1F("area_by_akt_jet", "Jet Areas When Found Individually for each akt jet", area_nBins, area_xMin, area_xMax);
    TH1F* areas_by_ca_jet = new TH1F("area_by_ca_jet", "Jet Areas When Found Individually for each ca reclustered jet", area_nBins, area_xMin, area_xMax);

    /*
    int pt_nBins = 200;
    double pt_xMin = 0;
    double pt_xMax = 200;
    TH1F* pt_akt = new TH1F("pt_akt", "pT of akt jets", pt_nBins, pt_xMin, pt_xMax);
    TH1F* pt_ca = new TH1F("pt_ca", "pT of ca reclustered jets", pt_nBins, pt_xMin, pt_xMax);
    */

    if (debug) cout << "Creating my event and background trees" << endl;
    my_event_tree met(event_file_name);
    my_background_trees mbt(background_file_names);

    if (debug) cout << "Calling getEntries()" << endl;
    int n_events = met.GetEntries();
    int n_backgrounds = mbt.GetEntries();

    fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));
    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme,fastjet::Best);
    fastjet::JetDefinition jet_def_ca(fastjet::cambridge_algorithm, Rparam, fastjet::E_scheme,fastjet::Best);

    int counter = 1000;
    for (int iEvent = 27000; iEvent < n_events; iEvent++) {

        if (iEvent >= counter) {
            cout << iEvent << " events processed" << endl;
            counter += 1000;
        }

        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max);

        vector<fastjet::PseudoJet> bg_particles = mbt.get_particles(iEvent);

        move(bg_particles.begin(), bg_particles.end(), back_inserter(event_particles));

        //fastjet::ClusterSequenceArea ca_akt_event(event_particles, jet_def_akt, area_def);
        fastjet::ClusterSequence ca_akt_event(event_particles, jet_def_akt);


        vector<fastjet::PseudoJet> akt_jets = ca_akt_event.inclusive_jets();
        if (debug) cout << "Looping through jets" << endl;
        for (fastjet::PseudoJet jet : akt_jets) {
            if (abs(jet.eta()) > jet_eta_max) continue;

            // For comparing area
            areas_by_event->Fill(jet.area(), met.bin_weight);

            fastjet::ClusterSequenceArea cs_akt_jet(jet.constituents(), jet_def_akt, area_def);

            areas_by_akt_jet->Fill(fastjet::sorted_by_pt(cs_akt_jet.inclusive_jets())[0].area(), met.bin_weight);

            fastjet::ClusterSequenceArea cs_ca_jet(jet.constituents(), jet_def_ca, area_def);
            areas_by_ca_jet->Fill(fastjet::sorted_by_pt(cs_ca_jet.inclusive_jets())[0].area(), met.bin_weight);
            
            // For comparing pT of akt jets and ca reclustered jets
            //pt_akt->Fill(jet.pt(), met.bin_weight);
            //pt_ca->Fill(fastjet::sorted_by_pt(fastjet::ClusterSequence(jet.constituents(), jet_def_ca).inclusive_jets())[0].pt(), met.bin_weight);

        } // end jet loop
        

    } // end event loop



    TCanvas *c1 = new TCanvas ("c1", "Canvas", 800, 600);

    
    c1->cd();
    areas_by_akt_jet->Draw("HIST");
    c1->Print(output_name + ".pdf(", "pdf");
    c1->Clear();

    c1->cd();
    areas_by_event->Draw("HIST");
    c1->Print(output_name + ".pdf", "pdf");
    c1->Clear();

    c1->cd();
    areas_by_ca_jet->Draw("HIST");
    c1->Print(output_name + ".pdf)", "pdf");
    c1->Clear();
    
    /*
    c1->cd();
    pt_akt->Draw("HIST");
    c1->Print("Comparing akt and ca pT.pdf(", "pdf");
    c1->Clear();

    c1->cd();
    pt_ca->Draw("HIST");
    c1->Print("Comparing akt and ca pT.pdf)", "pdf");
    c1->Clear();
    */


}