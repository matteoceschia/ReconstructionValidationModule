//! \file    ValidationModule.h
//! \brief   Example processing module for flreconstruct
//! \details Process a things object
#ifndef TESTMODULE_HH
#define TESTMODULE_HH
// Standard Library
// Third Party
//#include <boost/foreach.hpp>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TVector3.h"

// - Bayeux
#include "bayeux/dpp/base_module.h"
#include "bayeux/mctools/simulated_data.h"
#include "bayeux/genbb_help/primary_particle.h"
#include "bayeux/genbb_help/primary_event.h"
#include "bayeux/datatools/service_manager.h"
#include "bayeux/geomtools/manager.h"
#include "bayeux/geomtools/geometry_service.h"
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"
#include "bayeux/geomtools/geomtools.h"

// - Falaise
#include "falaise/snemo/datamodels/calibrated_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_data.h"
#include "falaise/snemo/datamodels/tracker_clustering_solution.h"
#include "falaise/snemo/datamodels/particle_track_data.h"


typedef struct ValidationEventStorage{
  // Quantities to histogram (h_)
  double h_total_calorimeter_energy_;
  double h_calo_energy_over_threshold_; // 50keV threshold to remove noise
  double h_calo_hit_time_separation_; // Between the first and last calo hits
  int h_calorimeter_hit_count_; // How many calorimeter hits?
  int h_calo_hits_over_threshold_; // How many calorimeter hits over threshold?
  int h_cluster_count_; // How many clusters with 3 or more hits?
  int h_track_count_; // How many reconstructed tracks?
  int h_negative_track_count_; // How many reconstructed tracks with negative curvature?
  int h_positive_track_count_; // How many reconstructed tracks with positive curvature?
  int h_associated_track_count_; // How many reconstructed tracks with an associated calorimeter?
  int h_geiger_hit_count_; // How many reconstructed tracker hits?
  double h_unassociated_calorimeter_energy_; // Summed calorimeter energy not associated to any track (in MeV)
  double h_unassociated_energy_over_threshold_; // Threshold is 50 keV
  double h_associated_calorimeter_energy_; // Summed calorimeter energy associated to any track (in MeV)
  double h_associated_energy_over_threshold_; // Threshold is 50 keV

  
  // For vector (v_) quantities you can have more than 1 entry per event
  // NOT implemented yet!
  std::vector<int> v_all_track_hit_counts_; // Vector of how many hits for ALL tracks (delayed or not)
  
  // For tracker maps: some values want to be summed over all events (t_), and some to be averaged (tm_)
  // All tm variables will need to be paired with a hit map, so that we can match the vector of values
  // to a vector of locations. This pairing must be defined in the config file.

  std::vector<int> t_cell_hit_count_; // map of cells that have been hit
  std::vector<double> tm_average_drift_radius_; // map of cells that have been hit
  // For calorimeter maps : some values want to be summed over all events (c_), and some to be averaged (cm_)
  // All cm variables will need to be paired with a hit map, so that we can match the vector of values
  // to a vector of locations. This pairing must be defined in the config file.
  std::vector<std::string> c_calorimeter_hit_map_;
  std::vector<double> cm_average_calorimeter_energy_;

}Validationeventstorage;


// This Project
class ValidationModule : public dpp::base_module {
  static const uint minHitsInCluster=3;
 public:
  //! Construct module
  ValidationModule();
  //! Destructor
  virtual ~ValidationModule();
  //! Configure the module
  virtual void initialize(const datatools::properties& myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);
  //! Process supplied data record
  virtual dpp::base_module::process_status process(datatools::things& workItem);
  //! Reset the module
  virtual void reset();
 private:
  TFile* hfile_;
  TTree* tree_;
  ValidationEventStorage validation_;

  // configurable data member
  std::string filename_output_;

  // geometry service
  const geomtools::manager* geometry_manager_; //!< The geometry manager

  void ResetVars();
  // You need to include these functions if you want to make detector maps
  int EncodeLocation(const snemo::datamodel::calibrated_tracker_hit & hit);
  std::string EncodeLocation(const snemo::datamodel::calibrated_calorimeter_hit & hit);

  // Macro which automatically creates the interface needed
  // to enable the module to be loaded at runtime
  DPP_MODULE_REGISTRATION_INTERFACE(ValidationModule);
};

#endif // TESTMODULE_HH


