#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftKiller.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"

#include "eventData.h"

#include "simpleTools.h"

using namespace std;

int main() {

    TString output_file_name = "akt softkilled jet particle phase space";

    TString event_file_name = "event_test.root";
    vector<TString> bg_file_names = {"backgrounds14ka.root", "backgrounds14kb.root"};

    // Jet stats
    double Rparam = .4;
    double part_eta_max = .9;
    double jet_eta_max = .9;

    vector<vector<double>> jet_pt_cuts = {
        {20, 120},
        {20, 40},
        {40, 60},
        {60, 80},
        {80, 100}, 
        {100, 120}

    };

    int r_numBins = 50;
    double r_min = 0;
    double r_max = Rparam + .05;

    int pt_numBins = 50;
    double pt_min = 0;
    double pt_max = 6;

    vector<vector<TH2F*>> all_hists;

    TString s = " ";
    for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
        vector<TH2F*> hists;
        hists.push_back(new TH2F("pythia_parts" + s + iCut, "R_{jet axis} and kt of pythia jet particles:" + s +jet_pt_cuts[iCut][0] + " < Jet p_{T} < " + jet_pt_cuts[iCut][1] + "; R_{jet_axis} 1/R; particle p_{T}", r_numBins, r_min, r_max, pt_numBins, pt_min, pt_max));
        hists.push_back(new TH2F("back_parts" + s + iCut, "R_{jet axis} and kt of background jet particles:" + s +jet_pt_cuts[iCut][0] + " < Jet p_{T} < " + jet_pt_cuts[iCut][1] , r_numBins, r_min, r_max, pt_numBins, pt_min, pt_max));

        all_hists.push_back(hists);

    }


    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best);
    fastjet::Selector jet_eta_selector = fastjet::SelectorAbsRapMax(jet_eta_max);
    fastjet::contrib::SoftKiller soft_killer(part_eta_max, .5);


    my_event_tree met(event_file_name);
    my_background_trees mbt(bg_file_names);


    int numEvents = met.GetEntries();
    int numBackgrounds = mbt.GetEntries();

    int counter = 1000;

    for (int iEvent = 0; iEvent < numEvents; iEvent++) {

        if (iEvent >= counter) {
            cout << iEvent << " events proccessed" << endl;
            counter += 1000;
        }


        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max);

        for (fastjet::PseudoJet &p : event_particles) { // index 0 for event particle, 1 for background particle
            p.set_user_index(0);
        }

        vector<fastjet::PseudoJet> bg_particles = mbt.get_particles(iEvent);
        for (fastjet::PseudoJet &p : bg_particles) { 
            p.set_user_index(1);
        }


        move (bg_particles.begin(), bg_particles.end(), back_inserter(event_particles));
        
        double pt_thresh = 0;
        vector<fastjet::PseudoJet> soft_killed_event;
        soft_killer.apply(event_particles, soft_killed_event, pt_thresh);
        event_particles = soft_killed_event;

        fastjet::ClusterSequence clust_seq(event_particles, jet_def_akt);
        vector<fastjet::PseudoJet> jets = jet_eta_selector(clust_seq.inclusive_jets());

        for (fastjet::PseudoJet jet : jets) {

            for (fastjet::PseudoJet p : jet.constituents()) {

                // find i cut based on jet pt
                for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {

                    if ((jet.pt() > jet_pt_cuts[iCut][0]) && (jet.pt() < jet_pt_cuts[iCut][1])) { // Find the right cuts

                        // p.user_index = 0 for pythia
                        //  particle, 1 for background particle
                        double delta = p.delta_R(jet);
                        all_hists[iCut][p.user_index()]->Fill(delta, p.pt(), met.bin_weight * (1 / delta)); // Need to weight by 1/delta because the amount of particles at a given radius increases by delta.
                    }
                }

            }

        } // end jet loop

    } // end event loop


    // print stuff out

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    TString title = "anti-kt jet particle phase space";
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + event_file_name;
    TString line3 = "Backgrounds from: " + print_vec(bg_file_names);

    TString line4 = "Jet radius: " + to_string_round(Rparam);
    TString line5 = "Particle eta max: " + to_string_round(part_eta_max);
    TString line6 = "Jet eta max: " + to_string_round(jet_eta_max);

    vector<TString> lines = {title, line1, line2, line3, line4, line5, line6};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_file_name + ".pdf(","pdf");
    c1->Clear();

    for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {

        all_hists[iCut][0]->Draw("COLZ");
        c1->Print(output_file_name + ".pdf", "pdf");
        c1->Clear();

        all_hists[iCut][1]->Draw("COLZ");
        if (iCut == jet_pt_cuts.size() - 1 ) c1->Print(output_file_name + ".pdf)", "pdf");
        else c1->Print(output_file_name + ".pdf", "pdf");
        c1->Clear();


    }

}