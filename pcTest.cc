//#include "fastjet/UserInfo.hh"
#include "pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"

#include "eventData.h"

using namespace Pythia8;

// From ecosia
void AddVerticalLine(TH1F* hist, double xValue) {
    TLine* line = new TLine(xValue, 0, xValue, hist->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2); // for dashed line
    line->Draw("SAME");
}

// From ecosia
void AddTextNote(TH1F* hist, double xValue, const char* text, Color_t color = kBlue) {
    // Add text near the x-axis
    TLatex* note = new TLatex(xValue, hist->GetMinimum() * 0.8, text);
    note->SetTextColor(color);
    note->SetTextSize(0.03);
    note->Draw("SAME");
}

int main() {

    TString eventFileName = "events1000x9.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";

    TString output_file_name = "jet vs pc tracks, all jets";

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

    double cone_area = M_PI * pow(Rparam, 2);

    //============================================================

    // Declaring graphs


    if(debug) cout << "Declaring histograms" << endl;
    // Specifically, tracks referring to background particles only
    int pt_numBins = 60;
    double pt_xMin = -.5;
    double pt_xMax = 20;
    double pt_step = (pt_xMax - pt_xMin) / pt_numBins;
    TH1F* jet_tracks_pt = new TH1F("jet_tracks_pt", "pT of Background Particles in Jet Cone; Track pT; frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    TH1F* pc_tracks_pt = new TH1F("pc_tracks_pt", "pT of Background Particles in PC Cone; Track pT; frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    
    TH1F* jet_jet_tracks_pt = new TH1F("jet_jet_tracks_pt", "pT of Background Particles in Jet Cone; Track pT; frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    TH1F* jet_pc_tracks_pt = new TH1F("jet_pc_tracks_pt", "pT of Background Particles in PC Cone; Track pT; frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    

    int area_numBins = 50;
    double area_xMin = 0;
    double area_xMax = 3;
    TH1F* akt_jet_area = new TH1F("akt_jet_area", "Area of akt jets; area; counts", area_numBins, area_xMin, area_xMax);
    // assuming area of pc is piR^2
    //TH1F* pc_area = new TH1F("pc_area", "Area of pc", area_numBins, area_xMin, area_xMax);

    TGraph* track_pt_ratio_uncorrected = new TGraph();
    track_pt_ratio_uncorrected->SetTitle("pt track-by-track ratio, jet/pc");
    TGraph* track_pt_ratio_area_corrected = new TGraph();

    //==========================================================

    // Prepare fastjet analysis
    fastjet::JetDefinition jet_def_akt(fastjet::antikt_algorithm, Rparam, fastjet::E_scheme, fastjet::Best);
    fastjet::AreaDefinition area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max));
    fastjet::Selector jet_eta_selector = fastjet::SelectorEtaRange(-1 * jet_eta_max, jet_eta_max);
    fastjet::Selector r_selector = fastjet::SelectorCircle(Rparam);

    //=============================================================

    //Begin event loop
    int numEvents = met.GetEntries();
    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }

    int n_jets = 0;
    int counter = 1000;

    

    for (int iEvent = 8500; iEvent < numEvents; iEvent++) {

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
            n_jets ++;
            jet_jet_tracks_pt->Reset("ICESM");
            jet_pc_tracks_pt->Reset("ICESM");

            akt_jet_area->Fill(jet.area());
            
            for (fastjet::PseudoJet &particle : jet.constituents()) {
                if (particle.user_index()==1) {
                    jet_tracks_pt->Fill(particle.pt());
                    jet_jet_tracks_pt->Fill(particle.pt());
                }
            }
            
            fastjet::PseudoJet pc_axis;
            pc_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));

            r_selector.set_reference(pc_axis);
            
            vector<fastjet::PseudoJet> pc_particles = r_selector(event_particles);
            for (fastjet::PseudoJet &particle : pc_particles) {
                if (particle.user_index() == 1) {
                    pc_tracks_pt->Fill(particle.pt());   
                    jet_pc_tracks_pt->Fill(particle.pt());   
                }
            }


            for (int iBin = 1; iBin <= pt_numBins; iBin++) {
                // Bin by bin, correct for the different area of the jet cone and the pc
                track_pt_ratio_area_corrected->AddPoint(pt_step*iBin, (jet.area() / cone_area) * (jet_tracks_pt->GetBinContent(iBin) / pc_tracks_pt->GetBinContent(iBin)));
            }


            if (leading_jet_only) break;
            
        } // end jet loop

        

    } // end event loop
    

    for (int iBin = 1; iBin <= pt_numBins; iBin++) {
        track_pt_ratio_uncorrected->AddPoint(pt_step*iBin, jet_tracks_pt->GetBinContent(iBin) / pc_tracks_pt->GetBinContent(iBin));
    }


    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    

    c1->cd();
    
    akt_jet_area->Draw("HIST");
    // Add pi R^2 line to akt_jet_area hist
    TLine* line = new TLine(cone_area, 0, cone_area, akt_jet_area->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2); // for dashed line
    line->Draw("SAME");

    c1->Print(output_file_name + ".pdf(", "pdf");
    c1->Clear();


    TLegend* legend = new TLegend(500, 700, 600, 800);
    legend->AddEntry(track_pt_ratio_uncorrected, "Uncorrected", "l");
    legend->AddEntry(track_pt_ratio_area_corrected, "Area corrected", "l");
    c1->cd();
    track_pt_ratio_uncorrected->Draw();
    track_pt_ratio_area_corrected->Draw("SAME");
    c1->Print(output_file_name + ".pdf", "pdf");
    c1->Clear();


    // When doing raw counts, pc_tracks_pt often has more because it has a larger area.
    /*
    vector<TH1F*> hists_to_normalize = {jet_tracks_pt, pc_tracks_pt};
    for (TH1F* hist : hists_to_normalize) {
        hist->Scale(1 / n_jets);
    }
    */

    vector<TH1F*> hists_to_print = {jet_tracks_pt, pc_tracks_pt};

    for (int iHist = 0; iHist < hists_to_print.size(); iHist++) {
        c1->cd();
        hists_to_print[iHist]->Draw("HIST");
        if (iHist == hists_to_print.size() - 1) c1->Print(output_file_name + ".pdf)", "pdf");
        else c1->Print(output_file_name + ".pdf", "pdf");
        c1->Clear();
    }

    // Now, we can go back through all the events, and look at what correction factor would account for upward background fluctuation recieved by the leading jet
    
    

}
