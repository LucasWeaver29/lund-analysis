#ifndef STORELUNDNEW_H
#define STORELUNDNEW_H

//#include "backSubTools.h"
#include "unifiedSubtractors"


struct lund_kin_vars {
    double delta, kt, jet_pt;
};

struct lund_coords {
    double log_inv_R, log_kt, jet_pt, bin_weight;
    int iBin;
};


// Making primary Lund plane from 
// https://hal.science/hal-01851158/document
class LundGroomer {

    public:
    
    // General
    double Rparam;
    bool leading_jet_only;
    bool debug = false;

    fastjet::JetDefinition jetDef_akt; // for jet identification
    //fastjet::JetDefinition jetDef_kt_recluster(fastjet::kt_algorithm, Rparam, recombScheme, strat); //  ???
    fastjet::Selector jet_eta_selector;

    EventSubtractor event_subtractor;
    JetSubtractor jet_subtractor;

    // For area based jet pT subtraction
    fastjet::AreaDefinition area_def;

    // Softdrop - Though SD (mine, all) can be imposed in lundFromFile
    double z_cut = .2;
    double beta = 0; // beta = 0 does just a normal z cut

    // Jet by jet constituent subtraction
    fastjet::GridMedianBackgroundEstimator bge; // GridMedianBackgroundEstimator is faster than JetMedianBackgroundEstimator and "performs equally well in nearly all cases"


    // RecursiveSoftDrop - fastjet contrib
    fastjet::contrib::RecursiveSoftDrop rsd;
    fastjet::contrib::SoftDrop sd;

    // So either can be used
    std::unique_ptr<fastjet::ClusterSequence> cs_ptr;
    std::unique_ptr<fastjet::ClusterSequenceArea> csa_ptr;
    
    double rho, jet_area;

    // constructor
    LundGroomer(double Rparam_in, double jet_eta_max_in, double part_eta_max_in, bool leading_jet_only_in = false):
    jetDef_akt(fastjet::antikt_algorithm, Rparam_in, fastjet::E_scheme, fastjet::Best),
    event_subtractor(part_eta_max_in),
    jet_subtractor(Rparam_in, part_eta_max_in),
    area_def(fastjet::active_area, fastjet::GhostedAreaSpec(part_eta_max_in)),
    bge(part_eta_max_in, .5), // particle eta max, grid spacing
    rsd(beta, z_cut, 5), // beta, z_cut, N (# of iterations)
    sd(beta, z_cut)
        {
        Rparam = Rparam_in;

        jet_eta_selector = fastjet::SelectorAbsRapMax(jet_eta_max_in);

        
        }

    
    
    // Takes event particles, int to count number of jets, int to record iCut, bg_sub_options.
    // Finds jets with akt. Reclusters those jets using CA.
    // returns a vector of the delta of the kinematic variables of those jets' splittings for the lund plane.
    vector<lund_kin_vars> get_kin_vars(vector<fastjet::PseudoJet> particles, unified_sub_options ops) {

        vector<lund_kin_vars> all_kinematic_vars;

        if (ops.event_sub != event_subtraction::null) {
            event_subtractor.subtract(&particles, ops.event_sub);
        }

        if (particles.empty()) return all_vars;


        if (debug) cout << "sotreLundNew LundGroomer: declustering particles" << endl;
        vector<fastjet::PseudoJet> jets;
        if (ops.jet_sub == jet_subtraction::Area_pT) { // If need area subtraction, need to record rho and cluster with area.
            bge.set_particles(particles);
            rho = bge.estimate().rho();
            csa_ptr = std::make_unique<fastjet::ClusterSequenceArea>(particles, jetDef_akt, area_def);
            jets = jet_eta_selector(csa_ptr->inclusive_jets());
        }
        else {
            cs_ptr = std::make_unique<fastjet::ClusterSequence>(particles, jetDef_akt);
            jets = jet_eta_selector(cs_ptr->inclusive_jets());
        }



        if (debug) cout << "Looping through akt jets" << endl;
        for (fastjet::PseudoJet& jet : jets) {  

            jet_counter ++;
            if (debug) {
                cout << "Jet num " << jet_counter << endl;
            }

            if (ops.jet_sub == jet_subtraction::Area_pT) jet_area = jet.area(); // Must get area before reclustering

            jet_subtractor.subtract_recluster(&jet, ops.jet_sub, particles);

            double jet_pt = jet.pt(); // Need to take pt after any subtractions have been applied. If no subtractions are applied, pt should be the same.
            if (sub_ops.jet_sub == jet_subtraction::Area_pT) {
                jet_pt -= jet_area * rho;
            }


            if (ops.weighted_jet_pt_hist != nullptr) {
                if(debug) cout << "Filling weighted jet pt hist: Jet pt =  " << jet.pt() << ", bin weight = " << ops.bin_weight << endl;
                ops.weighted_jet_pt_hist->Fill(jet_pt, ops.bin_weight);
            }

            
            if (debug) cout << "Beginning grooming + declustering" << endl;
            if (ops.groom_op == groom_option::SD_contrib_first) {
                jet = sd(jet);
            }

            fastjet::PseudoJet parent1, parent2;
            bool first = true;
            // For each jet, iteratively compare branchings
            while (jet.has_parents(parent1, parent2)) { 
                if (debug) cout << "New declustering branch" << endl;

                // In each case, identify the higher pt parent 
                fastjet::PseudoJet harder  = (parent1.pt() > parent2.pt()) ? parent1 : parent2;
                fastjet::PseudoJet softer  = (parent1.pt() > parent2.pt()) ? parent2 : parent1;
                
                double delta = harder.delta_R(softer);
                double kt = sin(delta) * softer.pt();
                
                if (first && ops.groom_op == groom_option::SD_mine_first) {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }


                if (ops.groom_op == groom_option::SD_mine_all) {
                    first = false;
                    // If softdrop condition is not met, discard lower pt subjet. Otherwise return kinematic variables and continue
                    if ((softer.pt() / (harder.pt() + softer.pt())) < z_cut * pow(delta/Rparam, beta)) {// soft drop condition
                        jet = harder;
                        continue;
                    }
                }

                lund_kin_vars vars = {.delta = delta, .kt = kt, .jet_pt = jet_pt};

                all_kinematic_vars.push_back(vars);

                jet = harder;
                first = false;
            } // end while loop

            if (leading_jet_only) break;

        } // end jet loop 
        return all_kinematic_vars;
    } // end decluster() def
};


#endif