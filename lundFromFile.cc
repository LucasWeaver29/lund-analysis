
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TCanvas.h"

#include "iostream"

#include <vector>

int main() {

    //std::vector<TString> lund_file_names = {"Pythia Lund, 2.3mil.root"};
    TString input_file = "Pythia Lund, 2.3mil.root";

    TString output_folder = "LundPlanes/2.3mil";
    TString output_file = "Pythia, from lundFromFile";

    bool debug = false;

    int jet_pt_min = 20;
    int jet_pt_max = 120;

    // For jet pt cuts of the lund plane
    std::vector<std::vector<double>> jet_pt_cuts = {
        {20, 40},
        {40, 60},
        {60, 80},
        {80, 100},
        {100, 120}
    };
    
    // For labeling lund plane
    TString lund_xAxis =  "ln((R0/delta)";
    TString lund_yAxis = "ln(kt)";
    
    double lund_xMin = 0;
    double lund_xMax = 5;
    int lund_nBinsX = 50;
    double lund_yMin = -4;
    double lund_yMax = 3.5;
    int lund_nBinsY = 50;

    TString space = " ";


    //========================================================
    // Setting up files and TTrees

    std::vector<TTree*> lund_coords_trees;
    std::vector<TH1F*> jet_pt_hists;

    //for (TString file_name : lund_file_names) {

        if(debug) std::cout << "Opening " << input_file << " as a ROOT file" << std::endl;

        TFile file(input_file, "Read");

        if (file.IsZombie()) {
            std::cout << "ERROR: Could not open ROOT file " << input_file << std::endl;
            return 1;
        }

        TTree* lund_coords_tree = (TTree*)file.Get("lund_coords_tree");
        if (!lund_coords_tree)  {
            std::cout << "ERROR: Could not find lund tree in Root file " << input_file << std::endl;
            return 1;
        }

        // Variables for tree branches
        double log_inv_R, log_kt, jet_pt, bin_weight;
        int iBin;

        if(debug) std::cout << "Setting branch addresses for lund_coords_tree" << std::endl;
        lund_coords_tree->SetBranchAddress("log_inv_R", &log_inv_R);
        lund_coords_tree->SetBranchAddress("log_kt", &log_kt);
        lund_coords_tree->SetBranchAddress("jet_pt", &jet_pt);
        lund_coords_tree->SetBranchAddress("bin_weight", &bin_weight);
        lund_coords_tree->SetBranchAddress("iBin", &iBin);



        //lund_coords_trees.push_back(lund_coords_tree);

        TTree* jet_pt_tree = (TTree*)file.Get("jet_pt_tree");
        if (!jet_pt_tree) {
            std::cout << "ERROR: Could not file jet pt tree in Root file " << input_file << std::endl;
            return 1;
        }

        TH1F* weighted_jet_pt_hist = nullptr;
        if(debug) std::cout << "Setting branch address for jet_pt_tree, and getting entry 0" << std::endl;
        jet_pt_tree->SetBranchAddress("jet_pt_hist", &weighted_jet_pt_hist);
        if(debug) std::cout << "Jet_pt_tree has " << jet_pt_tree->GetEntries() << " entries" << std::endl;
        jet_pt_tree->GetEntry(0);
        if(debug) std::cout << "Successfully got jet_pt_tree entry 0" << std::endl;
        

    //}


    //for (int iFile = 0; iFile < lund_file_names.size(); iFile++) {
        if(debug) std::cout << "lund_coords_tree has " << lund_coords_tree->GetEntries() << " entries" << std::endl;
        int num_points = lund_coords_tree->GetEntries();

        TH2F* inclusive_lund = new TH2F("inclusive_lund", "Inclusive Lund Plane;" + lund_xAxis + ";" + lund_yAxis, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax);

        std::vector<TH2F*> lund_cuts;
        if(debug) std::cout << "Adding TH2Fs to lund_cuts" << std::endl;
        for (int iCut = 0; iCut < jet_pt_cuts.size(); iCut++) {
            lund_cuts.push_back(new TH2F("lund_cut" + space + iCut, "Lund Jet p_{T} cut:" + space + jet_pt_cuts[iCut][0] + "< Jet p_{T} < " + jet_pt_cuts[iCut][1] + ";" + lund_xAxis + ";" + lund_yAxis, lund_nBinsX, lund_xMin, lund_xMax, lund_nBinsY, lund_yMin, lund_yMax));
        }


        for (int iPoint = 0; iPoint < num_points; iPoint++) {

            lund_coords_tree->GetEntry(iPoint);

            if ((jet_pt < jet_pt_min) || (jet_pt > jet_pt_max)) continue;

            inclusive_lund->Fill(log_inv_R, log_kt, bin_weight);

            // Find jet pt cut
            int iCut = 0;
            for (iCut; iCut < jet_pt_cuts.size(); iCut++) {
                if (jet_pt > jet_pt_cuts[iCut][0] && jet_pt < jet_pt_cuts[iCut][1]) break;
            }

            lund_cuts[iCut]->Fill(log_inv_R, log_kt, bin_weight);

        }


        // Normalize lund planes
        double nJets_weighted = 0;
        for (int iBin = weighted_jet_pt_hist->FindBin(jet_pt_min); iBin <= weighted_jet_pt_hist->FindBin(jet_pt_max); iBin++) {
            nJets_weighted += weighted_jet_pt_hist->GetBinContent(iBin);
        }


        double area = ((lund_xMax - lund_xMin)/lund_nBinsX) * ((lund_yMax - lund_yMin)/lund_nBinsY);
        double factor = 1/(area * nJets_weighted);
        inclusive_lund->Scale(factor);
        for (TH2F* lund_cut : lund_cuts) {
            lund_cut->Scale(factor);
        }
        
        //=========================================================
        // Print everything out
        //=========================================================


        TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

        // Adding a text information sheet with stats about this generation
        TString title = "Lund plane points from " + input_file;
        TString line1 = "Number of entries: " + std::to_string(num_points);
        TString line2 = "Jet p_{T} min = " + std::to_string(jet_pt_min);
        TString line3 = "Jet p_{T} max = " + std::to_string(jet_pt_max);
        TString line4 = "See " + input_file + "-Data_card.pdf for more info";



        std::vector<TString> lines = {title, line1, line2, line3, line4};
        c1->cd();
        for (int iLine = 0; iLine < lines.size(); ++iLine) {
            TLatex *text = new TLatex(.1, .8 - iLine * .05, lines[iLine]);
            text->SetTextSize(.04);
            text->Draw();
        }
        
        c1->Print(output_folder + "/"  + output_file + ".pdf(","pdf");
        c1->Clear();


        c1->cd();
        inclusive_lund->Draw("COLZ");
        c1->Print(output_folder + "/"  + output_file + ".pdf","pdf");
        c1->Clear();

        for (int iCut = 0; iCut < lund_cuts.size(); iCut++) {
            c1->cd();
            lund_cuts[iCut]->Draw("COLZ");
            (iCut == lund_cuts.size() - 1) ? 
                c1->Print(output_folder + "/"  + output_file + ".pdf)","pdf") :
                c1->Print(output_folder + "/"  + output_file + ".pdf","pdf");
            c1->Clear();
        }

    //}





}