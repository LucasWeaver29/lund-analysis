#include "pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"


// ROOT - files
//#include <TFile.h>
//#include <TTree.h>
//#include <TBranch.h>

#include <vector>
#include <iostream>
#include <filesystem>

// ROOT - histogram
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "eventData.h"
#include "backSubTools.h"

using namespace Pythia8;




struct lund_option {
    TString name;
    bg_sub_options sub_options;
};

int main() {
    
    //TString output_file_name = "9k, groomed";
    TString output_folder = "LundPlanes";
    TString subfolder = "debug weight issue";

    cout << "Make sure you have created the needed subfolder!" << endl;

    TString notes = "";
    //"z_cut = .2, beta = 3";


    vector<lund_option> pythia_lund_options = 
        {
            lund_option{.name = "From density Large", .sub_options = bg_sub_options{}},
            //lund_option{"Pythia, SD (Contrib's, first)", bg_sub_options{.groom_options = {"SD_contrib_first"}}},
            //lund_option{"Pythia, SD (mine, first)", bg_sub_options{.groom_options = {"SD_mine_first"}}},
            //lund_option{"Pythia, SD (mine, all)", bg_sub_options{.groom_options = {"SD_mine_all"}}},
        };

    vector<lund_option> bg_lund_options = 
        {
            //lund_option{.name = "Constituent Subtraction, PC particles", .sub_options = bg_sub_options{.pc_particles = true, .constit_sub = true}},
            //lund_option{"Background", bg_sub_options{}},
            //lund_option{"Background, SD (mine, all)", bg_sub_options{.groom_options = {"SD_mine_all"}}},
            //lund_option{"Constituent subtraction", bg_sub_options{.event_sub = "ConSub"}},
            //lund_option{"Constituent subtraction, SD (mine, all)", bg_sub_options{.event_sub = "ConSub", .groom_options = {"SD_mine_all"}}},
            //lund_option{"Perpendicular Cone (geometric)", bg_sub_options{.jet_sub = "MyPCGM"}},
            //lund_option{"Perpendicular Cone (geometric), SD (mine, all)", bg_sub_options{.jet_sub = "MyPCGM", .groom_options = {"SD_mine_all"}}}
            /*
            lund_option{"RSD (mine)", bg_sub_options{.jet_sub = "RSD_mine"}},
            lund_option{"RSD (mine), SD (mine, all)", bg_sub_options{.jet_sub = "RSD_mine", .groom_options = {"SD_mine_all"}}},
            lund_option{"RSD (contrib's)", bg_sub_options{.jet_sub = "RSD_contrib"}},
            lund_option{"RSD (mine), SD (mine, all)", bg_sub_options{.jet_sub = "RSD_mine", .groom_options = {"SD_mine_all"}}},
            */
        };


    TString eventFileName = "event_fromDensityLarge.root";
    //"event_scaled3Large.root";
    TString backgroundFileName = "thermalBackgrounds9000.root";
    //"thermalBackgrounds9000,etaMax=2.root";

    double ptMin = 20; // Minimum jet pT. 20 for Alice
    double ptMax = 120; // Maximum jet pt. 120 for Alice
    double part_eta_max = .9;
    double jet_eta_max = .5;
    double part_pt_min = .15;
    double Rparam = .4;
    bool leading_jet_only = false;
    if (part_eta_max - Rparam != jet_eta_max) cout << "WARNING: Rparam, part_eta_max, and jet_eta_max are not in agreement" << endl;

    bool debug = false;

    //==================================================================
    // Preparing to read from ROOT files
    my_event_tree met(eventFileName);
    my_background_tree mbt(backgroundFileName);
    // Checking if output subfolder exists
    //std::filesystem::path current_path = std::filesystem::current_path();
    //std::filesystem::path full_path = current_path / output_folder / subfolder

    //=======================================================================
    // Graphs 'n stuff
    // Declare 2D ROOT histograms
    TString xAxisName = "ln((R0/delta)";
    TString yAxisName = "ln(kt)";
    TString space = " ";
    double lund_xMin = 0;
    double lund_xMax = 5;
    int lund_nBinsX = 50;
    double lund_yMin = -4;
    double lund_yMax = 3.5;
    int lund_nBinsY = 50;

    int pt_nBins = 150;

    // jet pt cuts
    //vector<double> ptCutMins = {10, 20, 30, 50, 70, 100, 130, 160, 190};
    //vector<double> ptCutMaxs = {20, 30, 50, 70, 100, 130, 160, 190, 250};
    vector<double> ptCutMins = {20, 30, 50, 70, 90, 110, 130, 150};
    vector<double> ptCutMaxs = {30, 50, 70, 90, 110, 130, 150, 250};


    vector<vector<TH2F*>> pythia_lund_cuts;
    vector<TH1F*> pt_pythia;
    if (debug) cout << "Abt to assemble pythia_lund_cuts histograms" << endl;
    for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {

        vector<TH2F*> cuts;

        // Push back inclusive lund plane to index 0 of cuts
        cuts.push_back(new TH2F(pythia_lund_options[iOption].name, pythia_lund_options[iOption].name + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));

        // Add a histogram for tracking weighted jet pt counts to pt_bg
        pt_pythia.push_back(new TH1F("pt_pythia" + space + to_string(iOption), pythia_lund_options[iOption].name + ": jet pt spectrum; jet pt; weighted counts", pt_nBins, -.5, 299.5));
        pythia_lund_options[iOption].sub_options.jet_pt_hist = pt_pythia[iOption];

        
        // Add lund plane jet pt cuts
        for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            TString name = "pythia_lund_cuts" + pythia_lund_options[iOption].name + to_string(iCut);
            TString title = space + ptCutMins[iCut] + " < pT_{jet} < " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName;
            cuts.push_back(new TH2F(name, title, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }
        pythia_lund_cuts.push_back(cuts);

    }


    vector<TH2F*> bin_cuts;
    int numBins = met.GetNumBins();
    for (int iBin = 0; iBin < numBins; iBin++) {
        met.GetBinEntry(iBin);
        bin_cuts.push_back(new TH2F((TString)"binCut" + iBin, "pT Hat Bin:" + space + (*met.pt_hat_bin)[0] + " to " + (*met.pt_hat_bin)[1] + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
    
    }

    TH1F* bin_contributions_raw = new TH1F("bin_contributions_raw", "Raw contributions to Lund plane by pT hard bin; pT hard bin; N_{jets} raw", met.GetNumBins(), 0, met.GetNumBins());
    TH1F* bin_contributions_weighted = new TH1F("bin_contributions_weighted", "Weighted contributions to Lund plane by pT hard bin; pT hard bin; N_{jets} weighted", met.GetNumBins(), 0, met.GetNumBins());
    TH1F* bin_jets_raw = new TH1F("bin_jets_raw", "N_{jets} raw per pT hard bin; pT hard bin; N_{jets} raw", met.GetNumBins(), 0, met.GetNumBins());
    TH1F* bin_jets_weighted = new TH1F("bin_jets_weighted", "N_{jets} weighted per pT hard bin; pT hard bin; N_{jets} weighted", met.GetNumBins(), 0, met.GetNumBins());



    //vector<TString> bg_file_names;
    vector<vector<TH2F*>> bg_lund_cuts; // index 0 is inclusive
    vector<TH1F*> pt_bg;
    if (debug) cout << "Abt to assemble bg_lund_cuts histograms" << endl;
    for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {

        vector<TH2F*> cuts;

        // Push back inclusive lund plane to index 0 of cuts
        cuts.push_back(new TH2F(bg_lund_options[iOption].name, bg_lund_options[iOption].name + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));

        // Add a histogram for tracking weighted jet pt counts to pt_bg
        pt_bg.push_back(new TH1F("pt_bg" + space + to_string(iOption), bg_lund_options[iOption].name + ": jet pt spectrum; jet pt; weighted counts", pt_nBins, -.5, 299.5));
        bg_lund_options[iOption].sub_options.jet_pt_hist = pt_bg[iOption];

        // Add lund plane jet pt cuts
        for (int iCut = 0; iCut < ptCutMins.size(); ++iCut) {
            cuts.push_back(new TH2F("blc" + bg_lund_options[iOption].name + to_string(iCut), space + ptCutMins[iCut] + " < pT_{jet} < " + ptCutMaxs[iCut] + "; " + xAxisName + "; " + yAxisName, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }
        bg_lund_cuts.push_back(cuts);

    }
    
    // Set up fastjet analysis
    LundGroomer lund_groomer(ptMin, ptMax, leading_jet_only, Rparam, jet_eta_max, part_eta_max);
    lund_groomer.set_cuts(ptCutMins, ptCutMaxs);
    
    // =======================================================================

    // Reading events from ROOT
    int numEvents = met.GetEntries();
    //std::cout << "Number of events: " << numEvents << std::endl;

    /*
    int numBackgrounds = mbt.GetEntries();
    if (numBackgrounds < numEvents) {
        cout << "Warning: There are " << numEvents << " provided events but only" << numBackgrounds << " provided backgrounds." << endl;
        cout << "Lund plane creation will stop once backgrounds run out" << endl;
        numEvents = numBackgrounds;
    }
    */

    int counter = 1000;

    // Loop over events in ROOT TTree
    for (int iEvent = 0; iEvent < numEvents; ++iEvent) {
        

        if (debug) cout << "iEvent: " << iEvent << endl;

        if (iEvent > counter) {
            cout << counter << " events analyzed from ROOT file" << endl;
            counter += 1000;
        }

        if (debug) cout << "Calling met.get_particles" << endl;
        vector<fastjet::PseudoJet> event_particles = met.get_particles(iEvent, part_eta_max, part_pt_min); // iEvent, part_eta_max
        
        
        if (debug) cout << "Looping through pythia_lund_options" << endl;
        
        
        for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {
            int nJet = -1;
            if (debug) cout << "iOption " << iOption << endl;
            pythia_lund_options[iOption].sub_options.bin_weight = met.bin_weight;
            vector<lund_kin_vars> all_kin_vars = lund_groomer.get_kin_vars(event_particles, pythia_lund_options[iOption].sub_options);
            for (lund_kin_vars vars : all_kin_vars) {
                double logR0_R = log(Rparam/vars.delta);
                double logkt = log(vars.kt);
                pythia_lund_cuts[iOption][0]->Fill(logR0_R,logkt, met.bin_weight); // The inclusive pt lund
                pythia_lund_cuts[iOption][vars.num_cut + 1]->Fill(logR0_R,logkt, met.bin_weight); // the jet pt cut lund (which is at iCut + 1 in pythia_lund_cuts);
                bin_contributions_raw->Fill(met.event_iBin);
                bin_contributions_weighted->Fill(met.event_iBin, met.bin_weight);
                if (vars.num_jet > nJet) { // A new jet for this event
                    nJet = vars.num_jet;
                    bin_jets_raw->Fill(met.event_iBin);
                    bin_jets_weighted->Fill(met.event_iBin, met.bin_weight);
                }
                bin_cuts[met.event_iBin]->Fill(logR0_R,logkt, met.bin_weight);
            }
        }
        

        // Embedding in the background
        //vector<fastjet::PseudoJet> background_prtcls = mbt.get_particles(iEvent, part_pt_min);
        //move(background_prtcls.begin(), background_prtcls.end(), back_inserter(event_particles));

        if (debug) cout << "Looping through background_lund_options" << endl;
        for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {
            if (debug) cout << "iOption " << iOption << endl;
            bg_lund_options[iOption].sub_options.bin_weight = met.bin_weight;
            //bg_lund_options[iOption].sub_options.jet_pt_hist
            vector<lund_kin_vars> all_kin_vars = lund_groomer.get_kin_vars(event_particles, bg_lund_options[iOption].sub_options);
            for (lund_kin_vars vars : all_kin_vars) {
                double logR0_R = log(Rparam/vars.delta);
                double logkt = log(vars.kt);
                bg_lund_cuts[iOption][0]->Fill(logR0_R,logkt, met.bin_weight); // The inclusive pt lund
                bg_lund_cuts[iOption][vars.num_cut + 1]->Fill(logR0_R,logkt, met.bin_weight); // the jet pt cut lund (which is at iCut + 1 in pythia_lund_cuts);
            }
        }

    }  // End reading events from ROOT
    

// =========================================================
    // Save these plots in a root file, so they can be modified/resized without having to run the whole thing again
    /*
    TFile output_root_file(output_folder_name + root_output_name + ".root", "RECREATE");

    TH2F* lund_pointer = lund_pythia;

    TTree lund_tree("lund_tree", "Lund Tree");
    lund_tree.Branch("lund", &lund_pointer);
 
    // First entry is inclusive lund
    lund_tree.Fill();
    // subsequent entries are cuts
    for (int iLund = 0; iLund < lunds.size(); ++iLund) {
        lund_pointer = lunds[iLund];
        lund_tree.Fill();
    }

    lund_tree.Write();
    output_root_file.Close();
    */
//========================================================
    // Output files

    double area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY); // For normalization

    // Weighted analysis
    vector<TH1F*> weight_analysis_hists = {bin_contributions_raw, bin_contributions_weighted, bin_jets_raw, bin_jets_weighted};

    // Pythia ================================================
    for (int iOption = 0; iOption < pythia_lund_options.size(); ++iOption) {
        
        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

        // Adding a text information sheet with stats about this generation
        TString title = "Groomed Lund planes:; " + subfolder + "; " + pythia_lund_options[iOption].name;
        TString line1 = "Number of events: " + to_string(numEvents);
        TString line2 = "Events from: " + eventFileName;
        TString line3 = "Backgrounds from: " + backgroundFileName;
        TString line4 = "Leading Jet Only: " + bool2Str(leading_jet_only);

        TString line5 = "Jet pT min: " + to_string((int)ptMin);
        TString line6 = "Jet pT max: " + to_string((int)ptMax);

        TString line7 = "Jet radius: " + to_string_round(Rparam);
        TString line8 = "Particle eta max: " + to_string_round(part_eta_max);
        TString line9 = "Jet eta max: " + to_string_round(jet_eta_max);

        vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf(","pdf");
        c1->Clear();

        // Normalizing: Divide each bin by bin area and number of jets
        //cout << "nJets: " << nJets_bg[iFile] << endl;
        double nbgJets_weighted = 0;
        for(int iBin = pt_pythia[iOption]->FindBin(ptMin); iBin < pt_pythia[iOption]->FindBin(ptMax); ++iBin) {
            nbgJets_weighted += pt_pythia[iOption]->GetBinContent(iBin);
        }
        // Normalizing: Divide each bin by bin area and number of jets
        //area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY);
        double factor = 1/(area * nbgJets_weighted);
        for (TH2F* lund : pythia_lund_cuts[iOption]) {
            lund->Scale(factor);
        }
         
        for (int iLund = 0; iLund < pythia_lund_cuts[iOption].size(); ++iLund) {
            
            if (pythia_lund_cuts[iOption][iLund]->GetEntries() == 0) continue;

            c1->cd();
            pythia_lund_cuts[iOption][iLund]->Draw("COLZ");
            c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }

        for (TH1F* hist : weight_analysis_hists) {
            c1->cd();
            hist->Draw("HIST");
            c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }

        for (TH2F* lund : bin_cuts) {
            c1->cd();
            lund->Draw("COLZ");
            c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }
        

        c1->cd();
        c1->SetLogy(1);
        pt_pythia[iOption]->Draw("HIST");
        c1->Print(output_folder + "/" + subfolder + "/"  + pythia_lund_options[iOption].name + ".pdf)","pdf");
        c1->Clear();
        c1->SetLogy(0);

    } // End pythia file loop
    


    // Idealy, each bin is making the same number of raw contributions to the lund plane, because there is some number of contributions it takes to smooth out the plane
    // Multiplying the current number of events in each bin by the new multiplier will bring it to a number of events that gets desired_num_contributions to the lund plane
    int desired_num_contributions = 1;
    cout << "{";
    for (int iBin = 0; iBin < met.GetNumBins(); iBin++) {
        met.GetBinEntry(iBin);
        //print_bin(*met.pt_hat_bin);
        //cout << "Multiplier: " << desired_num_contributions / bin_contributions_raw->GetBinContent(iBin+1) << endl;
        cout << desired_num_contributions / bin_contributions_raw->GetBinContent(iBin+1) << ", ";
    }
    cout << "}" << endl;

    /*
    // bg Lunds ======================================================
    for (int iOption = 0; iOption < bg_lund_options.size(); ++iOption) {
        
        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);


        // Adding a text information sheet with stats about this generation
        TString title = "Groomed Lund planes:; " + subfolder + "; " + bg_lund_options[iOption].name;
        TString line1 = "Number of events: " + to_string(numEvents);
        TString line2 = "Events from: " + eventFileName;
        TString line3 = "Backgrounds from: " + backgroundFileName;
        TString line4 = "Leading Jet Only: " + to_string(leading_jet_only);

        TString line5 = "Jet pT min: " + to_string(ptMin);
        TString line6 = "Jet pT max: " + to_string(ptMax);

        TString line7 = "Jet radius: " + to_string(Rparam);
        TString line8 = "Particle eta max: " + to_string(part_eta_max);
        TString line9 = "Jet eta max: " + to_string(jet_eta_max);

        vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf(","pdf");
        c1->Clear();

        // Normalizing: Divide each bin by bin area and number of jets
        //cout << "nJets: " << nJets_bg[iFile] << endl;
        double nbgJets_weighted = 0;
        for(int iBin = pt_bg[iOption]->FindBin(ptMin); iBin < pt_bg[iOption]->FindBin(ptMax); ++iBin) {
            nbgJets_weighted += pt_bg[iOption]->GetBinContent(iBin);
        }
        // Normalizing: Divide each bin by bin area and number of jets
        //area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax-lund_yMin)/lund_nBinsY);
        double factor = 1/(area * nbgJets_weighted);
        for (TH2F* lund : bg_lund_cuts[iOption]) {
            lund->Scale(factor);
        }
         
        for (int iLund = 0; iLund < bg_lund_cuts[iOption].size(); ++iLund) {
            
            c1->cd();
            bg_lund_cuts[iOption][iLund]->Draw("COLZ");
            c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf","pdf");
            c1->Clear();
        }

        c1->cd();
        c1->SetLogy(1);
        pt_bg[iOption]->Draw("HIST");
        c1->Print(output_folder + "/" + subfolder + "/"  + bg_lund_options[iOption].name + ".pdf)","pdf");
        c1->Clear();
        c1->SetLogy(0);

    } // End file loop
    */

    // =======================================
    // Messing with it manually to directly subtract two lund planes
    /*
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    TString manual_file_name = "Constituent Subtraction, Jet - PC";

    // Adding a text information sheet with stats about this generation
    TString title = "Groomed Lund planes: " + subfolder + "; " + manual_file_name;
    TString line1 = "Number of events: " + to_string(numEvents);
    TString line2 = "Events from: " + eventFileName;
    TString line3 = "Backgrounds from: " + backgroundFileName;
    TString line4 = "Leading Jet Only: " + to_string(leading_jet_only);

    TString line5 = "Jet pT min: " + to_string(ptMin);
    TString line6 = "Jet pT max: " + to_string(ptMax);

    TString line7 = "Jet radius: " + to_string(Rparam);
    TString line8 = "Particle eta max: " + to_string(part_eta_max);
    TString line9 = "Jet eta max: " + to_string(jet_eta_max);

    vector<TString> lines = {title, notes, line1, line2, line3, line4, line5, line6, line7, line8, line9};

    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }

    c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf(","pdf");
    c1->Clear();


    for (int iLund = 0; iLund < bg_lund_cuts[0].size(); ++iLund) {
        c1->cd();
        TH2F* sub = (TH2F*)bg_lund_cuts[1][iLund]->Clone("sub");
        sub->Add(bg_lund_cuts[0][iLund], -1);
        sub->Draw("COLZ");
        if (iLund == bg_lund_cuts[0].size()-1) c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf)","pdf");
        else c1->Print(output_folder + "/" + subfolder + "/"  + manual_file_name + ".pdf","pdf");
        c1->Clear();
    }
    */
    // end manually messing


    // Since ROOT histograms point to the data they store, once the eventRootFile
    // closes they no longer have access to their data
    met.close_file();
    mbt.close_file();

    return 0;
}



