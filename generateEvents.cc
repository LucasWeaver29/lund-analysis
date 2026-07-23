// Generates events and stores them in ROOT

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "vector"

#include "TGraph.h"
#include "TCanvas.h"

using namespace Pythia8;

//===================================================================
// Various Bin Functions

// Used "refining bin function" data and dezmos
double bin_function(double x) {
    x = x - 2.5; // Since this was made to be left handed, but we actually want it at the center
    if (x < 65) return 2.4 * pow(10,9) * pow(x,-4.18) + 457;
    else return -.221 * pow(x, 1.51) + 538;
}
/*
double bin_function(double x) {
    return 7 * pow(10,9) * pow(x, -3.57) + 6000;
}
*/
/*
struct pt_hat_bin {
    int pt_hat_min, pt_hat_max, numEvents;
};
*/

// Bins according to bin function
vector<vector<int>> make_bins(int n_bins, int bin_pt_min, int bin_pt_max) {

    int event_multiplier = 2;


    vector<vector<int>> bins;
    double step = ((double)(bin_pt_max - bin_pt_min))/(n_bins);

    for (int iBin = 0; iBin < n_bins; iBin++) {
        int n_events = event_multiplier * (int) ((20.0 / n_bins ) * bin_function(bin_pt_min + (iBin + .5)*step)); // Using n_bins / 20 becuase scaled3 lite had 8 bins, and the bin equation is based off of it. it's iBin + .5 so that we're making the requst at the center of the bin
        vector<int> bin = {(int)(bin_pt_min + iBin*step), (int)(bin_pt_min + (iBin+1)*step), n_events};
        bins.push_back(bin);
    }

    return bins;
}


// Bins in inncrements of one. 
vector<vector<int>> make_small_bins(int bin_pt_min, int bin_pt_max) {

    vector<vector<int>> bins;
    for (int i = 0; i < bin_pt_max - bin_pt_min; i++) {
        vector<int> bin = {bin_pt_min + i, bin_pt_min + i + 1, 1000};
        bins.push_back(bin);
    }
    return bins;
}


void print_bins (vector<vector<int>> bins) {

    int totEvents = 0;

    for (int iBin = 0; iBin < bins.size(); iBin++) {
        cout << "Bin " << iBin << ": pt of " << bins[iBin][0] << " to " << bins[iBin][1] << ", " << bins[iBin][2] << " events." << endl;
        totEvents += bins[iBin][2];
    }

    cout << "Total # of events: " << totEvents << endl;
}

// Uses the multipliers required to ensure each pt hat bin has an equal number of contributions to the lund plane
vector<vector<int>> bins_from_density(int tot_num_events, int bin_width, double pt_hat_min, double pt_hat_max){
    // From groomed lund, raw number of jet contributions to lund plane per pt hat bin
    // have to match with make_small_bins(0,250);
    vector<double> multipliers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.166667, 0, 0.111111, 0.0833333, 0.0666667, 0.0526316, 0.0232558, 0.0178571, 0.012987, 0.00943396, 0.00806452, 0.00704225, 0.00440529, 0.00442478, 0.00408163, 0.00324675, 0.00275482, 0.00269542, 0.00218341, 0.00223214, 0.0018018, 0.00184843, 0.00153374, 0.00175131, 0.00123457, 0.00110375, 0.00118203, 0.00121065, 0.000948767, 0.000931099, 0.000862813, 0.00100604, 0.000803859, 0.000780031, 0.000704225, 0.000841751, 0.000726744, 0.000675219, 0.000668896, 0.000635324, 0.000646831, 0.000557414, 0.000627353, 0.000582751, 0.00056243, 0.000515198, 0.000509944, 0.000517331, 0.000489716, 0.00046729, 0.000474608, 0.000481232, 0.000426076, 0.000424628, 0.000453309, 0.0004363, 0.00044843, 0.00044603, 0.00041425, 0.000421585, 0.000392927, 0.000370508, 0.000366569, 0.000385951, 0.000381971, 0.000360881, 0.000359712, 0.000363108, 0.000379075, 0.000358809, 0.000344353, 0.000384172, 0.000366838, 0.000366166, 0.000354484, 0.000354862, 0.000326797, 0.000324992, 0.000301114, 0.000344947, 0.000331016, 0.000319489, 0.000311526, 0.000305998, 0.000304321, 0.000295508, 0.000303398, 0.000283527, 0.000299222, 0.000275786, 0.000292141, 0.000293772, 0.000282486, 0.000274876, 0.000289017, 0.000281611, 0.000290782, 0.000283447, 0.000270197, 0.000277008, 0.000280269, 0.000277008, 0.000263783, 0.000256345, 0.00024919, 0.000271003, 0.000270709, 0.000261712, 0.000265887, 0.00026096, 0.000241196, 0.000261165, 0.000265393, 0.000262055, 0.000238265, 0.000244021, 0.000256739, 0.000258398, 0.000250627, 0.000252461, 0.000263019, 0.000252781, 0.000240038, 0.000243665, 0.000255624, 0.000241429, 0.000237699, 0.000248571, 0.000240385, 0.000255428, 0.000249004, 0.000241896, 0.000244738, 0.000252845, 0.000251004, 0.000228519, 0.000243487, 0.000247341, 0.000230097, 0.000236518, 0.000217014, 0.000238322, 0.000235571, 0.000237699, 0.000230097, 0.000240385, 0.000233427, 0.000238039, 0.000227894, 0.000236239, 0.000225124, 0.000238095, 0.000234577, 0.000225683, 0.000244081, 0.000266454, 0.000236239, 0.000247036, 0.000269324, 0.000237812, 0.000234962, 0.000236967, 0.000234632, 0.000227998, 0.000261097, 0.000242601, 0.000253293, 0.000248633, 0.000255624, 0.000249128, 0.000256739, 0.00026455, 0.000236183, 0.000244081, 0.000250815, 0.000248324, 0.000246245, 0.000260078, 0.000247463, 0.000270416, 0.000253807, 0.000263158, 0.000237925, 0.000251256, 0.000256213, 0.000263435, 0.000259808, 0.000251383, 0.000270783, 0.0002531, 0.000259538, 0.000264061, 0.000283768, 0.000257136, 0.000276932, 0.000304229, 0.000272628, 0.000276091, 0.000275406, 0.000273598, 0.000275103, 0.000265745, 0.000272851, 0.000258732, 0.000291121, 0.000260756, 0.000265534, 0.000295596, 0.000284414, 0.000280899, 0.000282167, 0.000289184, 0.000262812, 0.000280741, 0.00030525, 0.000319387, 0.000302663, 0.000309789, 0.000278707, 0.000319285, 0.000302939, 0.000274574, 0.000312989, 0.000308833, 0.000308928, 0.000282247, 0.00028777, 0.000325203, 0.000309885, 0.000322581, 0.000301296, 0.000292483, 0.000310945, 0.000316256, 0.000317259, 0.00031736, 0.00033036, 0.00032175};

    double event_density_factor = (tot_num_events / 21500.0) * 150000.0; // 21500 was how many events were made when 150000.0 was the factor, and we were going 20 < jet pt hat < 250

    vector<vector<int>> bins;    
    
    int width_counter = 0;
    double event_density = 0;

    for (int pt_hat = pt_hat_min; pt_hat < pt_hat_max; pt_hat++) {

        if (width_counter == bin_width) {
            int n_events = (int) (event_density * event_density_factor);
            vector<int> bin = {pt_hat - bin_width, pt_hat, n_events};
            bins.push_back(bin);

            width_counter = 0;
            event_density = 0;
        }
        
        event_density += multipliers[pt_hat];
        width_counter++;
        
    }

    return bins;


    // For bins of size 1
    /*
    vector<vector<int>> bins;
    int total_events = 0;
    // multipliers are in increments of one, as in multipliers[1] is from 1 < pt hat < 2.
    for (int pt_hat = 0; pt_hat < multipliers.size(); pt_hat++) {
        if (pt_hat < pt_hat_min) continue;
        if (pt_hat > pt_hat_max) break;
        
        int nEvents = (int) (150000.0 * multipliers[pt_hat]);
        vector<int> bin = {pt_hat, pt_hat + 1, nEvents};
        bins.push_back(bin);

        total_events+=nEvents;
    }

    cout << total_events << " total events to be generated." << endl;

    return bins;
    */
}


int main() {
    
    // Create a ROOT output file
    TString filename = "event_test.root";
    TFile file(filename, "RECREATE");

    // number of events per bin
    // do 5000

    double etaMax = 5;
    
    // Give bounds for each bin
    
    // using make_bins
    //std::vector<vector<int>> pt_hat_bins = bins_from_density(500000, 10, 20, 250); // total desired events, bin width, pt hat min, pt hat max
    //print_bins(pt_hat_bins);

    // For test
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 10000},
        {30, 50, 8000},
        {50, 70, 4000},
        {70, 100, 2000},
        {100, 250, 4000},
    };


    // Zoltan's binning
    /*
    std::vector<vector<int>> pt_hat_bins = {
            {20, 30, 1000000},
            {30, 50, 800000},
            {50, 70, 400000},
            {70, 100, 200000},
            {100, 130, 100000},
            {130, 180, 100000},
            {180, 250, 100000}
        };
    */

    // for small test, for making sure storeLund and lundFromFile are working
    /*
    std::vector<vector<int>> pt_hat_bins = {
        {20, 70, 100},
        {70, 120, 100}
    };
    */

    /*
    // Adding more statistic in jet pt > 70 range, for "event_moreHighPtJets.root"
    std::vector<vector<int>> pt_hat_bins = {
        {50, 70, 500000},
        {70, 90, 500000},
        {90, 110, 500000},
        {110, 130, 400000},
        {130, 150, 300000},
        {150, 250, 200000}
    };
    */

    // Scaled3 Lite
    /*
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 160000},
        {30, 50, 32500},
        {50, 70, 13500},
        {70, 90, 9500},
        {90, 110, 7750},
        {110, 130, 7000},
        {130, 150, 6000},
        {150, 250, 6000}
    };
    */
    
    
    /*
    // Scaled3 
    std::vector<vector<int>> pt_hat_bins = {
        {20, 30, 64000},
        {30, 50, 12000},
        {50, 70, 5100},
        {70, 90, 3800},
        {90, 110, 2800},
        {110, 130, 2700},
        {130, 150, 2500},
        {150, 250, 2500}
    };
    
    */


    /*
    vector<vector<int>> pt_hat_bins = {
        {10, 20, 1000},
        {20, 30, 1000},
        {30, 50, 1000},
        {50, 70, 1000},
        {70, 90, 1000},
        {90, 110, 1000},
        {110, 130, 1000},
        {130, 150, 1000},
        {150, 250, 1000}
    };
    */
    
    




    //vector<double> ptHatMins = {10, 20, 30, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {20, 30, 50, 70, 90, 110, 130, 150, 250};
    //vector<double> ptHatMins = {10, 50, 70, 90, 110, 130, 150};
    //vector<double> ptHatMaxs = {50, 70, 90, 110, 130, 150, 250};
    


    // Create a TTree
    TTree tree("events", "Pythia Events Tree");
    
    // Define variables to store event data
    int ptHatMin, ptHatMax; 
    int event_iBin;
    double pTHat;
    int numParticles;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> energy;
    std::vector<double> eta;
    std::vector<int> id;
    
    // Create branches for the tree
    tree.Branch("ptHatMin", &ptHatMin, "ptHatMin/I");
    tree.Branch("ptHatMax", &ptHatMax, "ptHatMax/I");
    tree.Branch("event_iBin", &event_iBin, "event_iBin/I");
    //tree.Branch("crossSection", &crossSection, "crossSection/D");
    tree.Branch("numParticles", &numParticles, "numParticles/I");
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("eta", &eta);
    tree.Branch("pdgId", &id);
    tree.Branch("pTHat", &pTHat);


    vector<int> pt_hat_bin;
    TTree ptHatBins_tree("ptHatBins", "Pythia pT Hat Bins Tree");
    ptHatBins_tree.Branch("pt_hat_bin", &pt_hat_bin);
    for (int iBin = 0; iBin < pt_hat_bins.size(); iBin++) {
        pt_hat_bin = pt_hat_bins[iBin];
        ptHatBins_tree.Fill();
    }
    ptHatBins_tree.Write();

    // For testing how pythia.info.sigmaGen() evolves over time. 
    TGraph* xsecs = new TGraph();
    xsecs->SetTitle("Pythia event cross sections over events");
    
    // To store bin weights such that they can be added back to the tree entries later
    vector<double> bin_weights(pt_hat_bins.size());

    // =================================================================

    // Initialize Pythia
    Pythia pythia;
    pythia.readString("Beams:eCM = 5360");
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numbershowEvent = 0");

    int tot_events = 0;

    // Bin loop
    for (int iBin = 0; iBin < pt_hat_bins.size(); ++iBin) {

        ptHatMin = pt_hat_bins[iBin][0];
        ptHatMax = pt_hat_bins[iBin][1];
        int nEvents = pt_hat_bins[iBin][2];
        tot_events++;

        pythia.settings.readString("PhaseSpace:pTHatMin = " + to_string(ptHatMin));
        pythia.settings.readString("PhaseSpace:pTHatMax = " + to_string(ptHatMax));

        // Proton - proton collision is default
        pythia.init();

        event_iBin = iBin;

        // Event loop
        for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
            
            if (!pythia.next()) continue; // Skip failed events
            
            // Fill event data
                    
            //crossSection = pythia.info.sigmaGen();
            //weight = pythia.info.sigmaGen() / nEvents;
            xsecs->AddPoint(tot_events, pythia.info.sigmaGen());

            // Clear vectors for new event
            px.clear();
            py.clear();
            pz.clear();
            energy.clear();
            eta.clear();
            id.clear();

            // Fill particle data
            for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
                
                // Only take final, charged particles in correct eta bin
                if (not pythia.event[iPart].isFinal()) continue;

                if (abs(pythia.event[iPart].eta()) > etaMax) continue;

                if (not pythia.event[iPart].isCharged()) continue;
                

                const Pythia8::Particle& particle = pythia.event[iPart];
                //std::cout << "    Particle " << iPart << ":" << std::endl;
                //std::cout << "      Pythia says px is " << pythia.event[iPart].px() << std::endl;
                px.push_back(particle.px());
                //std:cout << "      px vector says its "<< px.back() << endl;
                py.push_back(particle.py());
                pz.push_back(particle.pz());
                energy.push_back(particle.e());
                eta.push_back(particle.eta());
                id.push_back(particle.id());
                
                // Doesn't work bc iPart doesn't match px idexes
                /*
                std::cout << "      px: " << px[iPart] << std::endl;
                std::cout << "      py: " << py[iPart] << std::endl;
                std::cout << "      pz: " << pz[iPart] << std::endl;
                std::cout << "      energy: " << energy[iPart] << std::endl;
                std::cout << "      pdgId: " << id[iPart] << std::endl;
                */
            }

            numParticles = px.size();

            // Fill the tree
            tree.Fill();
        } // end event loop
        
        bin_weights[iBin] = pythia.info.sigmaGen() / nEvents;
        
    }


    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    c1->cd();
    c1->SetLogy(1);
    xsecs->Draw();
    c1->Print("Pythia Cross Sections Over Time.pdf","pdf");
    c1->Clear();


    // Write the tree to the file
    tree.Write();
    //weight_tree.Write();
    file.Close();




    // Reopen file to add bin_weight
    TFile *reopened_file = new TFile(filename, "UPDATE");
    if (reopened_file->IsZombie()) {
        cout << "Error: Could not reopen file " << filename << endl;
        return 1;
    }

    // Reaccess the tree
    TTree* reopened_tree = (TTree*)reopened_file->Get("events");
    if (!reopened_tree) {
        std::cout << "Error: Could not find 'events' tree in ROOT file" << std::endl;
    }

    // Add new branch for bin weights
    double bin_weight = 0.0;
    TBranch *weight_branch = reopened_tree->Branch("bin_weight", &bin_weight, "bin_weight/D");

    // Set branch address for iBin
    int reopened_iBin;
    reopened_tree->SetBranchAddress("event_iBin", &reopened_iBin);

    int nEntries = reopened_tree->GetEntries();
    for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        reopened_tree->GetEntry(iEntry);
        bin_weight = bin_weights[reopened_iBin];
        weight_branch->Fill();
    }

    // Save with overwrite flag
    reopened_tree->Write("", TObject::kOverwrite);
    reopened_file->Close();


    return 0;
}
