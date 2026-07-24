#ifndef UNIFIEDSUBTRACTORS_H
#define UNIFIEDSUBTRACTORS_H

#include "backSubTools.h"

using namespace std;



enum class event_subtraction {
    ConSub,
    ConSubRho,
    SoftKill,
    null
};


enum class jet_subtraction {
    // Requires different cluster sequence
    Area_pT,
    
    // Done with constituents
    MyPCGM,
    MyPCkT,
    PC_parts,
    // Requires first clusterings with CA, then getting constituents
    RSD_mine,

    // Done with an already CA reclustered jet
    Filter,
    Prune,
    RSD_contrib,
    null
};

enum class groom_option {
    SD_contrib_first,
    SD_mine_first,
    SD_mine_all,
    null
};

struct unified_sub_options {
    TString name;
    event_subtraction event_sub = event_subtraction::null;
    jet_subtraction jet_sub = jet_subtraction::null;
    groom_option groom_op = groom_option::null;
    TH1F* weighted_jet_pt_hist = nullptr;
    double bin_weight;
    bool embed_in_bg;
};

class EventSubtractor {
    public:

    double part_eta_max;

    bool debug = false;

    // Jet by jet constituent subtraction
    fastjet::GridMedianBackgroundEstimator bg_estimator; // GridMedianBackgroundEstimator is faster than JetMedianBackgroundEstimator and "performs equally well in nearly all cases"
    fastjet::contrib::ConstituentSubtractor constituent_subtractor;

    // SoftKiller
    fastjet::contrib::SoftKiller soft_killer;


    EventSubtractor(double part_eta_max_in):
    bg_estimator(part_eta_max_in, .5), // Particle eta max, grid spacing
    constituent_subtractor(&bg_estimator),
    soft_killer(part_eta_max_in, .5) // particle eta max, grid spacing
    {
        part_eta_max = part_eta_max_in;

        // Set up subtractor for jet by jet constituent subtraction
        constituent_subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
        constituent_subtractor.set_max_distance(0.3);  // R_max for ghost-particle pairing
        constituent_subtractor.set_max_eta(part_eta_max_in); // These two are for whole-event mode
        constituent_subtractor.initialize(); //
    }


    void subtract(vector<fastjet::PseudoJet>* particles, event_subtraction sub) {
        if (debug) cout << "Beginning event wide subtractions" << endl;
        
        if (sub == event_subtraction::ConSub) {
            bg_estimator.set_particles(*particles);
            //vector<fastjet::PseudoJet> temp_particles = constituent_subtractor.subtract_event(*particles);
            //*particles = temp_particles;
            *particles = constituent_subtractor.subtract_event(*particles);
            if (particles->empty() && debug) cout << "No remaining particles after contrib constituent subtraction" << endl;
        }
        else if (sub == event_subtraction::ConSubRho) {
            bg_estimator.set_particles(*particles);
            //double rho = rho_estimator.rho();
            fastjet::contrib::ConstituentSubtractor manual_subtractor(bg_estimator.rho());
            *particles = manual_subtractor.subtract_event(*particles);
            if (particles->empty() && debug) cout << "No remaining particles after bigger rho constituent subtraction" << endl;
        }
        else if (sub == event_subtraction::SoftKill) {
            vector<fastjet::PseudoJet> soft_killed_event;
            double pt_thresh = 0; // pt threshold for killed particles
            soft_killer.apply(*particles, soft_killed_event, pt_thresh);
            *particles = soft_killed_event;
        }

    } // end subtract_event() def
        

};




class JetSubtractor {

    public:

    double Rparam;

    bool debug = false;

    fastjet::JetDefinition jet_def_recluster;

    // Softdrop, rsd, etc
    double z_cut = .2;
    double beta = 0;

    // my pc
    my_pc_subtractor pc_subtractor;

    //filter 
    double Rfilt = .3;
    double nfilt = 3;
    fastjet::Filter filter;

    // Pruning
    double prune_zcut = .2;
    double Rcut_factor = .2;
    fastjet::Pruner pruner;

    // RecursiveSoftDrop - fastjet contrib
    fastjet::contrib::RecursiveSoftDrop rsd;
    int rsd_n_iterations = 5;

    // Need to keep the cluster sequence alive after using it
    std::unique_ptr<fastjet::ClusterSequence> cs_ptr;
    
    // If neccessary for area based subtraction - Actually, this needs to be done with the whole event at once, and should be done in the principal script
    //std::unique_ptr<fastjet::ClusterSequenceArea> csa_ptr;

    //fastjet::AreaDefinition area_def;


    // constructor
    JetSubtractor(double Rparam_in, double part_eta_max_in):
    jet_def_recluster(fastjet::cambridge_algorithm, Rparam_in + .4, fastjet::E_scheme, fastjet::Best),
    pc_subtractor(Rparam_in, part_eta_max_in),
    filter(Rfilt, fastjet::SelectorNHardest(nfilt)),
    pruner(fastjet::cambridge_algorithm, prune_zcut, Rcut_factor),
    rsd(beta, z_cut, rsd_n_iterations)
    //area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max_in))
    {
        Rparam = Rparam_in;
    }



    // For when jets need to be reclustered with ca
    void subtract_recluster(fastjet::PseudoJet* jet, jet_subtraction sub, vector<fastjet::PseudoJet> particles) {

        vector<fastjet::PseudoJet> constituents;

        if (sub == jet_subtraction::MyPCGM) {
            constituents = pc_subtractor.geometric_subtract_constit(*jet, particles);
        }
        else if (sub == jet_subtraction::MyPCkT) {
            constituents = pc_subtractor.kt_subtract_constit(*jet, particles);
        }
        else if (sub == jet_subtraction::PC_parts) {
            constituents = pc_subtractor.get_pc_particles(*jet, particles);
        }
        else if (sub == jet_subtraction::RSD_mine) {
            fastjet::ClusterSequence reclusterSeq_temp(jet->constituents(), jet_def_recluster);
            fastjet::PseudoJet reclustered_jet = sorted_by_pt(reclusterSeq_temp.inclusive_jets())[0];
            constituents = recursive_soft_drop_constit(reclustered_jet, Rparam);
        }
        else constituents = jet->constituents();


        if (constituents.empty()) {
            if (debug) cout << "No remainaing constituents after jet subtraction" << endl;
            fastjet::PseudoJet empty_jet(0,0,0,0);
            constituents.push_back(empty_jet);
        } 
        

        if (debug) cout << "JetSubtractor: reclustering with new jet constituents" << endl;

        /*
        if (sub == jet_subtraction::Area_pT) {
            csa_ptr = std::make_unique<fastjet::ClusterSequenceArea>(constituents, jet_def_recluster, area_def);
            //cs_ptr = nullptr;
            *jet = fastjet::sorted_by_pt(csa_ptr->inclusive_jets())[0];
        }
        */
        
        cs_ptr = std::make_unique<fastjet::ClusterSequence>(constituents, jet_def_recluster);
        *jet = fastjet::sorted_by_pt(cs_ptr->inclusive_jets())[0];
        
        
        

        // Next, the things that work only on an already CA reclustered jet, and return a new jet
        if  (sub == jet_subtraction::Filter) {
            *jet = filter(*jet);
        }
        else if (sub == jet_subtraction::Prune) {
            *jet = pruner(*jet);
        }
        else if (sub == jet_subtraction::RSD_contrib) {
            *jet = rsd(*jet);
        }

    }



};





#endif