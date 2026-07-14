#include <TFile.h>
#include <TTree.h>
#include <iostream>

// From ecosia
// Remeber that iBins will now be meaningless if the two event file were generated with different binning!
// This also will not copy the pt_hat_bins tree.
void appendTTrees(const char* file1_name, const char* file2_name, const char* output_file_name) {

    TFile* file1 = new TFile(file1_name, "READ");
    if (file1->IsZombie()) {
        std::cout << "Error: Could not open root file " << file1_name << std::endl;
    }

    TFile* file2 = new TFile(file2_name, "READ");
    if (file2->IsZombie()) {
        std::cout << "Error: Could not open root file " << file2_name << std::endl;
    }

    TFile* output_file = new TFile(output_file_name, "RECREATE");

    // Get the TTrees from the files
    TTree* tree1 = (TTree*)file1->Get("events");
    TTree* tree2 = (TTree*)file2->Get("events");

    TTree* mergedTree = tree1->CloneTree(0);
    std::cout << "Copied format of Tree1" << std::endl;

    mergedTree->CopyEntries(tree1);
    std::cout << "Copied entries of tree1 into mergedTree" << std::endl;
    mergedTree->CopyEntries(tree2);

    mergedTree->Write();
    output_file->Close();

    file1->Close();
    file2->Close();

}   


int main() {
    appendTTrees("event_2.3mil.root", "event_moreHighPtJets.root", "event_5mil.root");
}
