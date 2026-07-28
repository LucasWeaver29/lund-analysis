//#include "fastjet/UserInfo.hh"
//#include "pythia8/Pythia.h"

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

#include "simpleTools.h"

using namespace std;

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

    TString eventFileName = "event_test.root";
    vector<TString> backgroundFileNames = {"backgrounds14ka.root", "backgrounds14kb.root"};

    TString output_file_name = "jet vs pc tracks";

    bool debug = false;

    if (debug) cout << "Declaring my_event_tree" << endl;
    my_event_tree met(eventFileName);
    if (debug) cout << "Declaring my_background_tree" << endl;
    my_background_trees mbt(backgroundFileNames);


    double Rparam = .4;

    double part_eta_max = .9;
    double jet_eta_max = .5;

    double part_pt_min = .15;

    double jet_pt_min = 0;
    double jet_pt_max = 500;

    bool leading_jet_only = false;

    double cone_area = M_PI * pow(Rparam, 2);

    //============================================================

    // Declaring graphs


    if(debug) cout << "Declaring histograms" << endl;
    // Specifically, tracks referring to background particles only
    int pt_numBins = 60;
    double pt_xMin = 0;
    double pt_xMax = 3.5;
    double pt_step = (pt_xMax - pt_xMin) / pt_numBins;
    TH1F* jet_tracks_pt = new TH1F("jet_tracks_pt", "pT of Background Particles in Jet Cone; Track pT; #frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    TH1F* pc_tracks_pt = new TH1F("pc_tracks_pt", "pT of Background Particles in PC Cone; Track pT; #frac{1}{N_jets} dpT", pt_numBins, pt_xMin, pt_xMax);
    TH1F* pc_tracks_pt_area_corrected = new TH1F("pc_tracks_pt_area_corrected", "pT of Background Particles in PC, area corrected; Track pT; #frac{#Area_jet}{#Area_PC} #frac{1}{N_jets} d#p_T", pt_numBins, pt_xMin, pt_xMax);
    
    int area_numBins = 50;
    double area_xMin = 0;
    double area_xMax = .75;
    TH1F* akt_jet_area = new TH1F("akt_jet_area", "Area of akt jets; area; counts", area_numBins, area_xMin, area_xMax);
    // assuming area of pc is piR^2
    //TH1F* pc_area = new TH1F("pc_area", "Area of pc", area_numBins, area_xMin, area_xMax);

    /*
    TGraph* track_pt_ratio_uncorrected = new TGraph();
    track_pt_ratio_uncorrected->SetTitle("pt track-by-track ratio, pc/jet; track pt");
    track_pt_ratio_uncorrected->SetLineColor(kRed);

    TGraph* track_pt_ratio_area_corrected = new TGraph();
    track_pt_ratio_area_corrected->SetLineColor(kBlue);
    */

    TGraph* cfs = new TGraph();
    cfs->SetTitle("PC Correction Factors; track pt; multiplier");
    cfs->SetLineColor(kBlue);

    TGraph* cfs_no_area_correction = new TGraph();
    cfs_no_area_correction->SetLineColor(kRed);
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
        jets = sorted_by_pt(jets);
        
        if (debug) cout << "Starting jet loop" << endl;
        for (fastjet::PseudoJet &jet : jets) {
            n_jets ++;
            
            if (jet.pt() > jet_pt_max) continue;
            if (jet.pt() < jet_pt_min) break; // Break because jets are in descending pt order
            //temp_jet_tracks_pt->Reset("ICESM");
            //temp_pc_tracks_pt->Reset("ICESM");

            akt_jet_area->Fill(jet.area(), met.bin_weight);
            
            for (fastjet::PseudoJet &particle : jet.constituents()) {
                if (particle.user_index()==1) {
                    jet_tracks_pt->Fill(particle.pt(), met.bin_weight);
                    //temp_jet_tracks_pt->Fill(particle.pt());
                }
            }
            
            fastjet::PseudoJet pc1_axis, pc2_axis;
            pc1_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() + M_PI/2 , 2*M_PI));
            pc2_axis.reset_PtYPhiM(0,jet.eta(), fmod(jet.phi() - M_PI/2 , 2*M_PI));

            vector<fastjet::PseudoJet> pc_axes = {pc1_axis, pc2_axis};
            for (fastjet::PseudoJet axis : pc_axes) {
                
                r_selector.set_reference(axis);
            
                vector<fastjet::PseudoJet> pc_particles = r_selector(event_particles);
                for (fastjet::PseudoJet &particle : pc_particles) {
                    if (particle.user_index() == 1) { // Factors of .5 since we're doing two perpendicular cones
                        pc_tracks_pt->Fill(particle.pt(), met.bin_weight * .5);   
                        //temp_pc_tracks_pt->Fill(particle.pt(), .5);
                        pc_tracks_pt_area_corrected->Fill(particle.pt(), met.bin_weight * (jet.area() / cone_area) * .5);      
                    }
                }
            }
            

            if (leading_jet_only) break;
            
        } // end jet loop

        

    } // end event loop
    

    for (int iBin = 1; iBin <= pt_numBins; iBin++) {
        if ((jet_tracks_pt->GetBinContent(iBin) == 0) || (pc_tracks_pt->GetBinContent(iBin) == 0)) continue;
        cfs_no_area_correction->AddPoint(pt_step * iBin, jet_tracks_pt->GetBinContent(iBin) / pc_tracks_pt->GetBinContent(iBin));
        //track_pt_ratio_uncorrected->AddPoint(pt_step * iBin, pc_tracks_pt->GetBinContent(iBin) / jet_tracks_pt->GetBinContent(iBin));
    }

    TH1F* correction_factors_hist = new TH1F("correction_facotrs_hist", "PC Correction Factors per Track pT Bin; track pT; Correction Factor", pt_numBins, pt_xMin, pt_xMax);
    for (int iBin = 1; iBin <= pt_numBins; iBin++) {
        if (pc_tracks_pt_area_corrected->GetBinContent(iBin) == 0 || jet_tracks_pt->GetBinContent(iBin) == 0) continue;    
        //track_pt_ratio_area_corrected->AddPoint(pt_step * iBin, pc_tracks_pt_area_corrected->GetBinContent(iBin) / jet_tracks_pt->GetBinContent(iBin));
        cfs->AddPoint(pt_step * iBin, jet_tracks_pt->GetBinContent(iBin) / pc_tracks_pt_area_corrected->GetBinContent(iBin));
        correction_factors_hist->Fill(pt_step * iBin, jet_tracks_pt->GetBinContent(iBin) / pc_tracks_pt_area_corrected->GetBinContent(iBin));
    }

    // Store correction factor hist for later use
    TFile cf_output_file("PC_Correction_Factors.root", "RECREATE");
    TTree* pc_cf_tree = new TTree("pc_cf_tree", "PC Correction Factors Tree");
    pc_cf_tree->Branch("cf_hist", &correction_factors_hist);
    pc_cf_tree->Fill();
    pc_cf_tree->Write();
    cf_output_file.Close();



    //=====================================================================================================
    // Printing stuff...
    //=====================================================================================================


    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    
    TString title = "Perpendicular Cone Correction Factor Analysis";
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = "Backgrounds from: " + print_vec(backgroundFileNames);
    TString line4 = "Leading Jet Only: " + bool2Str(leading_jet_only);

    TString line5 = "Jet pT min: " + to_string((int)jet_pt_min);
    TString line6 = "Jet pT max: " + to_string((int)jet_pt_max);

    TString line7 = "Jet radius: " + to_string_round(Rparam);
    TString line8 = "Particle eta max: " + to_string_round(part_eta_max);
    TString line9 = "Jet eta max: " + to_string_round(jet_eta_max);

    vector<TString> lines = {title, line1, line2, line3, line4, line5, line6, line7, line8, line9};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_file_name + ".pdf(","pdf");
    c1->Clear();



    c1->cd();
    
    akt_jet_area->Draw("HIST");
    // Add pi R^2 line to akt_jet_area hist
    TLine* line = new TLine(cone_area, 0, cone_area, akt_jet_area->GetMaximum());
    line->SetLineColor(kRed);
    line->SetLineStyle(2); // 2 for dashed line
    line->Draw("SAME");

    c1->Print(output_file_name + ".pdf", "pdf");
    c1->Clear();



    // Need to find ymin and ymax of ratio TGraphs together so that both will be in view
    double graph_ymin = .6;
    double graph_ymax = 1.7;
    
    c1->cd();

    cfs->GetHistogram()->SetMaximum(graph_ymax);
    cfs->GetHistogram()->SetMinimum(graph_ymin);


    cfs->Draw();
    cfs_no_area_correction->Draw("SAME");

    TLegend* legend = new TLegend(.7, .7, .9, .9);
    legend->AddEntry(cfs, "Correction Factors", "l");
    legend->AddEntry(cfs_no_area_correction, "Correction Factors, no area correction", "l");
    legend->Draw();
    
    c1->Print(output_file_name + ".pdf", "pdf");
    c1->Clear();


    // When doing raw counts, pc_tracks_pt often has more because it has a larger area.
    /*
    vector<TH1F*> hists_to_normalize = {jet_tracks_pt, pc_tracks_pt};
    for (TH1F* hist : hists_to_normalize) {
        hist->Scale(1 / n_jets);
    }
    */

    vector<TH1F*> hists_to_print = {jet_tracks_pt, pc_tracks_pt, pc_tracks_pt_area_corrected};

    for (int iHist = 0; iHist < hists_to_print.size(); iHist++) {
        c1->cd();
        hists_to_print[iHist]->Draw("HIST");
        if (iHist == hists_to_print.size() - 1) c1->Print(output_file_name + ".pdf)", "pdf");
        else c1->Print(output_file_name + ".pdf", "pdf");
        c1->Clear();
    }

    // Now, we can go back through all the events, and look at what correction factor would account for upward background fluctuation recieved by the leading jet
    
    met.close_file();
    mbt.close_file();

}