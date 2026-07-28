
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"

#include "iostream"

#include <vector>

#include "simpleTools.h"

using namespace std;


int main() {

    TString overlay_output_name = "lund overlay, no grooming, best";


    //vector<TString> lund_file_names = {"Pythia Lund, 2.3mil.root"};
    vector<TString> input_files = {
        // Not groomed
        
        //"Pythia, from event_Zoltans",
        //"Embedded, from event_Zoltans",
        //"My PC, from event_Zoltans",
        //"RSD (mine), from event_Zoltans",
        //"ConSub, from event_Zoltans",
        //"Embedded, Area Sub, from event_Zoltans",
        //"SoftKill, from event_Zoltans"

        
        
        // GROOMED with SD (mine, all)
        /*
        "Pythia, SD (mine, all), from event_Zoltans",
        "RSD (mine), SD (mine, all), from event_Zoltans",
        "Embedded, SD (mine, all), from event_Zoltans",
        "SoftKill, SD (mine, all), from event_Zoltans",
        "My PC, SD (mine, all), from event_Zoltans",
        "ConSub, SD (mine, all), from event_Zoltans"
        */
        
        

        // From event_test, on my pc
        
        "Pythia, SD (mine, all), from event_test",
        "Embedded, SD (mine, all), from event_test",
        "Embedded, Area Sub, from event_test",
        "ConSub, SD (mine, all), from event_test", 
        "My PC Geometric Sub, SD (mine, all), from event_test"
        //"MyPC, with cf, from event_test"
        
    }; 
    // dont' include ".root", that's done automatically
    //"Pythia, SD (mine, all), from event_Zoltans"
    //TString input_files_folder = "LundPlanes/event_Zoltans/rootFiles/";
    // The first lund plane cut of the first file will be used as the reference

    vector<TString> file_short_names = {};
    for (TString name : input_files) {
        TString short_name = name;
        file_short_names.push_back(short_name.ReplaceAll(", from event_Zoltans", ""));
    }

    TString output_folder = "LundPlanes";

    bool debug = false;

    // For jet pt cuts of the lund plane
    vector<vector<double>> jet_pt_cuts = {
        {20, 120}, // first cuts is the inclusive plane
        {20, 40},
        {40, 60},
        {60, 80},
        {80, 100},
        {100, 120}
    };

    // cuts along log(kt) and log(R0/R)
    vector<double> kt_cuts = {-2, -1, 0, 2}; // Lines of ln(kt) = x along which to take cuts
    vector<double> delta_cuts = {.25, .5, 1, 1.5};
    int cut_width = 2; // The number of adjecent, parallel rows or columns to include in a cut

    
    // For labeling lund plane
    TString lund_xAxis =  "ln((R0/delta)";
    TString lund_yAxis = "ln(kt)";
    
    double lund_xMin = 0;
    double lund_xMax = 5;
    int lund_nBinsX = 50;
    double lund_yMin = -4;
    double lund_yMax = 3.5;
    int lund_nBinsY = 50;

    double lund_area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax - lund_yMin)/lund_nBinsY);


    TString space = " ";



    //========================================================
    // Setting up files and TTrees
    
    vector<vector<TH2F*>> all_lunds;


    for (int iFile = 0; iFile < input_files.size(); iFile++) {

        TString input_file = input_files[iFile];

        vector<TH2F*> lunds;
        
        if(debug) cout << "Opening " << input_file << " as a ROOT file" << endl;
        TFile file(input_file + ".root", "Read");

        if (file.IsZombie()) {
            cout << "ERROR: Could not open ROOT file " << input_file << endl;
            return 1;
        }

        TTree* lund_coords_tree = (TTree*)file.Get("lund_coords_tree");
        if (!lund_coords_tree)  {
            cout << "ERROR: Could not find lund tree in Root file " << input_file << endl;
            return 1;
        }

        // Variables for tree branches
        double log_inv_R, log_kt, jet_pt, bin_weight;
        int iBin;

        if(debug) cout << "Setting branch addresses for lund_coords_tree" << endl;
        lund_coords_tree->SetBranchAddress("log_inv_R", &log_inv_R);
        lund_coords_tree->SetBranchAddress("log_kt", &log_kt);
        lund_coords_tree->SetBranchAddress("jet_pt", &jet_pt);
        lund_coords_tree->SetBranchAddress("bin_weight", &bin_weight);
        lund_coords_tree->SetBranchAddress("iBin", &iBin);

        //===============================
        // Get the weighted jet pt histogram from the file's jet_pt_Tree
        TTree* jet_pt_tree = (TTree*)file.Get("jet_pt_tree");
        if (!jet_pt_tree) {
            cout << "ERROR: Could not find jet pt tree in Root file " << input_file << endl;
            return 1;
        }

        TH1F* weighted_jet_pt_hist = nullptr;
        if(debug) cout << "Setting branch address for jet_pt_tree, and getting entry 0" << endl;
        jet_pt_tree->SetBranchAddress("jet_pt_hist", &weighted_jet_pt_hist);
        jet_pt_tree->GetEntry(0);

        //lunds.weighted_jet_pt_hist = (TH1F*)weighted_jet_pt_hist->Clone();
        
        //===============================


        if(debug) cout << "lund_coords_tree has " << lund_coords_tree->GetEntries() << " entries" << endl;
        int num_points = lund_coords_tree->GetEntries(); // Points to graph on the lund plane

        
        if(debug) cout << "Adding TH2Fs to jet_pt_cuts" << endl;
        for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
            TH2F* l = new TH2F(input_file + "lund_cut" + iCut, jet_pt_cuts[iCut][0] + space + "< Jet p_{T} < " + jet_pt_cuts[iCut][1] + ";" + lund_xAxis + ";" + lund_yAxis, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax);
            l->SetDirectory(0);
            lunds.push_back(l);
        }

        if(debug) cout << "Adding points to lund planes" << endl;
        // Add the points in the lund_coords_tree to the lund planes
        for (int iPoint = 0; iPoint < num_points; iPoint++) {

            lund_coords_tree->GetEntry(iPoint);

            // Find jet pt cut
            for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
                if (jet_pt > jet_pt_cuts[iCut][0] && jet_pt < jet_pt_cuts[iCut][1]) {
                    lunds[iCut]->Fill(log_inv_R, log_kt, bin_weight); 
                }
            }


        }

        //===========================================================
        // Normalize lund planes
        if(debug) cout << "Normalizing lund planes" << endl;
        for (int iLund = 0; iLund < lunds.size(); iLund++) { 
            double nJets_weighted = 0;
            for (int iBin = weighted_jet_pt_hist->FindBin(jet_pt_cuts[iLund][0]); iBin <= weighted_jet_pt_hist->FindBin(jet_pt_cuts[iLund][1]); iBin++) {
                nJets_weighted += weighted_jet_pt_hist->GetBinContent(iBin);
            }
            lunds[iLund]->Scale(1/(lund_area * nJets_weighted));
        }
        
        all_lunds.push_back(lunds);

        //=========================================================
        // Print individual files out
        //=========================================================

        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

        // Adding a text information sheet with stats about this generation
        TString title = "Lund plane points from " + input_file;
        TString line1 = "Number of entries: " + to_string(num_points);
        TString line2 = "See " + input_file + "-Data_card.pdf for more info";


        vector<TString> lines = {title, line1, line2};
        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/"  + input_file + ".pdf(","pdf");
        c1->Clear();

        for (int iLund = 0; iLund < lunds.size(); iLund++) {
            c1->cd();
            lunds[iLund]->Draw("COLZ");
            //(iLund == lunds.size() - 1) ? 
                //c1->Print(output_folder + "/"  + input_file + ".pdf)","pdf") :
                c1->Print(output_folder + "/"  + input_file + ".pdf","pdf");
            c1->Clear();
        }

        c1->cd();
        c1->SetLogy(1);
        weighted_jet_pt_hist->Draw("HIST");
        c1->Print(output_folder + "/"  + input_file + ".pdf)","pdf");
        c1->Clear();

        file.Close();

    }


    //=========================================================
    // Do kt and R cuts to compare different lunds, save to lund overlay.pdf
    //=========================================================

    TH2F *reference_lund = (TH2F*)all_lunds[0][0]->Clone(); // First cut of first file is used as reference
    reference_lund->SetStats(false);
    reference_lund->SetTitle(";;");
    
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    // Adding a text information sheet with stats about this lund overlay
    TString file_list = "";
    for (int iFile = 0; iFile < input_files.size(); iFile++) {
        if (iFile == input_files.size() - 1) file_list += input_files[iFile] + ".";
        else file_list += input_files[iFile] + ", ";
    }

    TString title = "Overlay of lunds from " + file_list;
    
    
    vector<TString> lines = {title};
    c1->cd();
    for (int iLine = 0; iLine < lines.size(); ++iLine) {
        TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
        text->SetTextSize(.04);
        text->Draw();
    }
    
    c1->Print(output_folder + "/" + overlay_output_name + ".pdf(","pdf");
    c1->Clear();

    // Do kt and R cuts for each jet pt cut. i_JetPtCut = 0 corresponds to the inclusive lund plane, jet pt cuts are indexed from 1
    for (int i_pTCut = 0; i_pTCut < all_lunds[0].size(); i_pTCut++) { 
        if(debug) cout << "Beginning jet pt cut " << i_pTCut << endl;
        
        // kt cuts
        for (int i_kt = 0; i_kt < kt_cuts.size(); i_kt++) { 
            if(debug) cout << "Beginning kt cut " << i_kt << endl;

            TLegend *legend = new TLegend(.1,.7,.28,.9);

            c1->cd();

            // To ensure no histogram is cut off
            double display_yMax = 0;
            //double display_yMin = 1000; // a big number so TMath::Min(GetMax(), display_yMin) logic works
            TH1D* first_hist = nullptr; // The first hist drawn sets the display axes


            for (int iFile = 0; iFile < all_lunds.size(); iFile++) { // Iterate through the lunds from the different input files. These are the entries in the vector all_lunds

                int cut_iBin = all_lunds[iFile][i_pTCut]->GetYaxis()->FindBin(kt_cuts[i_kt]); // Find the bin corresponding to the kt cut

                if (debug) cout << "Projecting hist for kt_cut" << endl;
                TH1D* kt_cut = all_lunds[iFile][i_pTCut]->ProjectionX("kt_cut" + iFile, cut_iBin - cut_width, cut_iBin + cut_width, "e"); // Include cut_width bins on either side of the exact ln(kt) value
                kt_cut->Scale(1.0/(2*cut_width +1)); // to account for the multiple rows summed for each cut
                display_yMax = max(display_yMax, kt_cut->GetMaximum());
                //display_yMin = min(display_yMin, kt_cut->GetMinimum());

                if (iFile==0) {
                    kt_cut->SetTitle(all_lunds[iFile][i_pTCut]->GetTitle() + space + ": ln(kt) =" + space + kt_cuts[i_kt] + "; ln(R0/R); weighted counts");
                    first_hist = kt_cut;
                    kt_cut->Draw("HIST" "PLC");
                    kt_cut->SetStats(false); // no stat box

                }
                else kt_cut->Draw("HIST" "SAME" "PLC");
            
                legend->AddEntry(kt_cut, file_short_names[iFile], "l");

            }
            
            first_hist->GetYaxis()->SetRangeUser(0, display_yMax * 1.1);
            c1->Update();

            // Add the reference
            TPad* pad = new TPad("pad", "pad", .8, .75, .95, .95);
            pad->Draw();
            pad->cd(); // Switch to smaller pad
            TH2F* ref_copy = (TH2F*)all_lunds[0][i_pTCut]->Clone();
            ref_copy->SetStats(false);
            ref_copy->SetTitle(";;");
            ref_copy->Draw("COLZ A25");

            TLine* line = new TLine(lund_xMin, kt_cuts[i_kt], lund_xMax, kt_cuts[i_kt]);
            line->SetLineColor(kRed);
            line->SetLineWidth(2); 
            line->Draw("SAME");
            
            c1->cd(); // switch back to main canvas

            legend->Draw();
            c1->Print(output_folder + "/" + overlay_output_name + ".pdf", "pdf");
            c1->Clear();
        }

        // delta cuts
        for (int iDelta = 0; iDelta < delta_cuts.size(); iDelta++) { 
            if (debug) cout << "Beginning delta cut " << iDelta << endl;

            TLegend *legend = new TLegend(.1,.7,.28,.9);

            c1->cd();

            double display_yMax = 0;
            //double display_yMin = 1000; 
            TH1D* first_hist = nullptr;
            
            for (int iFile = 0; iFile < all_lunds.size(); iFile++) { // Iterate through the lunds from the different input files. These are the entries in the vector all_lunds

                int cut_iBin = all_lunds[iFile][i_pTCut]->GetXaxis()->FindBin(delta_cuts[iDelta]);
                
                TH1D* delta_cut = all_lunds[iFile][i_pTCut]->ProjectionY("detla_cut" + iFile, cut_iBin - cut_width, cut_iBin + cut_width, "e");
                delta_cut->Scale(1.0/(2*cut_width + 1));
                display_yMax = max(display_yMax, delta_cut->GetMaximum());
                //display_yMin = min(display_yMin, delta_cut->GetMinimum());


                if (iFile == 0) {
                    delta_cut->SetTitle(all_lunds[iFile][i_pTCut]->GetTitle() + space + ": ln(R0/R) =" + space + delta_cuts[iDelta] + "; ln(kt); weighted, normalized counts");
                    first_hist = delta_cut;
                    delta_cut->Draw("HIST" "PLC");
                    delta_cut->SetStats(false); // no stat box
                }
                else delta_cut->Draw("HIST" "SAME" "PLC");
            
                legend->AddEntry(delta_cut, file_short_names[iFile], "l");
            }
            
            first_hist->GetYaxis()->SetRangeUser(0, display_yMax * 1.05);
            c1->Update();

            // Add the reference
            TPad* pad = new TPad("pad", "pad", .8, .75, .95, .95);
            pad->Draw();
            pad->cd(); // Switch to smaller pad
            TH2F* ref_copy = (TH2F*)all_lunds[0][i_pTCut]->Clone();
            ref_copy->SetStats(false);
            ref_copy->SetTitle(";;");
            ref_copy->Draw("COLZ A25");

            TLine* line = new TLine(delta_cuts[iDelta], lund_yMin, delta_cuts[iDelta], lund_yMax);
            line->SetLineColor(kRed);
            line->SetLineWidth(2); // for dashed line
            line->Draw("SAME");
            
            c1->cd(); // switch back to main canvas

            legend->Draw();
            if ((iDelta == delta_cuts.size() - 1) && (i_pTCut == all_lunds[0].size() - 1)) c1->Print(output_folder + "/" + overlay_output_name + ".pdf)", "pdf");
            else c1->Print(output_folder + "/" + overlay_output_name + ".pdf", "pdf");
            c1->Clear();

        }

    }



}