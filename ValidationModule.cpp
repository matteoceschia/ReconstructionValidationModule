#include "ValidationModule.h"
const double LOW_ENERGY_LIMIT=0.050; // 50 keV

int mainWallHitType=1302;
int xWallHitType=1232;
int gammaVetoHitType=1252;

using namespace std;


DPP_MODULE_REGISTRATION_IMPLEMENT(ValidationModule,"ValidationModule");
ValidationModule::ValidationModule() : dpp::base_module()
{
  filename_output_="Validation.root";
}

ValidationModule::~ValidationModule() {
  if (is_initialized()) this->reset();
}

void ValidationModule::initialize(const datatools::properties& myConfig,
                                   datatools::service_manager& flServices,
                                   dpp::module_handle_dict_type& /*moduleDict*/){

  // Look for services
  if (flServices.has("geometry")) {
    const geomtools::geometry_service& GS = flServices.get<geomtools::geometry_service> ("geometry");
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }
  // Extract the filename_out key from the supplied config, if
  // the key exists. datatools::properties throws an exception if
  // the key isn't in the config, so catch this if thrown and don't do
  // anything
  try {
    myConfig.fetch("filename_out",this->filename_output_);
  } catch (std::logic_error& e) {
  }

  // Use the method of PTD2ROOT to create a root file with just the branches we need for the Validation analysis


  hfile_ = new TFile(filename_output_.c_str(),"RECREATE","Output file of Simulation data");
  hfile_->cd();
  tree_ = new TTree("Validation","Validation");
  tree_->SetDirectory(hfile_);
 
  // Some basic counts
  tree_->Branch("h_calorimeter_hit_count",&validation_.h_calorimeter_hit_count_);
  tree_->Branch("h_calo_hits_over_threshold",&validation_.h_calo_hits_over_threshold_);
  tree_->Branch("h_cluster_count",&validation_.h_cluster_count_);
  tree_->Branch("h_track_count",&validation_.h_track_count_);
  tree_->Branch("h_negative_track_count",&validation_.h_negative_track_count_);
  tree_->Branch("h_positive_track_count",&validation_.h_positive_track_count_);
  tree_->Branch("h_associated_track_count",&validation_.h_associated_track_count_);
  tree_->Branch("h_geiger_hit_count",&validation_.h_geiger_hit_count_);
  tree_->Branch("v_all_track_hit_counts",&validation_.v_all_track_hit_counts_);
  
  // Energies and calo times
  tree_->Branch("h_total_calorimeter_energy",&validation_.h_total_calorimeter_energy_);
  tree_->Branch("h_calo_energy_over_threshold",&validation_.h_calo_energy_over_threshold_);
  tree_->Branch("h_unassociated_calorimeter_energy",&validation_.h_unassociated_calorimeter_energy_);
  tree_->Branch("h_associated_calorimeter_energy",&validation_.h_associated_calorimeter_energy_);
  tree_->Branch("h_unassociated_energy_over_threshold",&validation_.h_unassociated_energy_over_threshold_);
  tree_->Branch("h_associated_energy_over_threshold",&validation_.h_associated_energy_over_threshold_);
  tree_->Branch("h_calo_hit_time_separation",&validation_.h_calo_hit_time_separation_);
  
  // Tracker maps
  tree_->Branch("t_cell_hit_count",&validation_.t_cell_hit_count_);
  tree_->Branch("tm_average_drift_radius",&validation_.tm_average_drift_radius_);

  // Calo maps
  tree_->Branch("c_calorimeter_hit_map",&validation_.c_calorimeter_hit_map_);
  tree_->Branch("c_average_calorimeter_energy",&validation_.cm_average_calorimeter_energy_);
  
  this->_set_initialized(true);
}
//! [ValidationModule::Process]
dpp::base_module::process_status
ValidationModule::process(datatools::things& workItem) {
  
  
  // declare internal variables to mimic the ntuple variables, names are same but in camel case
  double totalCalorimeterEnergy=0;
  double energyOverThreshold=0;
  double unassociatedEnergy=0;
  double unassocOverThreshold=0;
  double timeDelay=-1;
  int clusterCount=0;
  int trackCount=0;
  int alphaCount=0;
  int associatedTrackCount=0;
  int caloHitCount=0;
  int nCalHitsOverLowLimit=0;
  double earliestCaloHit=-1.;
  double latestCaloHit=-1.;
  int geigerHitCount=0;
  int negativeTrackCount=0;
  int positiveTrackCount=0;
  std::vector<int> allTrackHitCounts;


  // We need to run this before we start populating vectors. Put all your vectors in this function to clear them
  ResetVars();
  
  // Grab calibrated data bank
  // Calibrated data will only be present in reconstructed files,
  // so wrap in a try block



  try {
    const snemo::datamodel::calibrated_data& calData = workItem.get<snemo::datamodel::calibrated_data>("CD");

    
    if (calData.has_calibrated_calorimeter_hits())
      {
        const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & calHits=calData.calibrated_calorimeter_hits();
        for (snemo::datamodel::calibrated_data::calorimeter_hit_collection_type::const_iterator   iHit = calHits.begin(); iHit != calHits.end(); ++iHit) {
          const snemo::datamodel::calibrated_calorimeter_hit & calHit = iHit->get();
          
          // Write to the calorimeter map
          validation_.c_calorimeter_hit_map_.push_back(EncodeLocation(calHit));
          
          double energy=calHit.get_energy() ;
          totalCalorimeterEnergy += energy;
          if (energy > LOW_ENERGY_LIMIT)
          {
            energyOverThreshold +=energy;
            ++nCalHitsOverLowLimit;
          }
          double hitTime=calHit.get_time();
          if (hitTime > latestCaloHit)
          {
            latestCaloHit = hitTime;
          }
          if (hitTime < earliestCaloHit || earliestCaloHit <= 0 )
          {
            earliestCaloHit=hitTime;
          }
          ++caloHitCount;
        }
      }
    
    timeDelay=latestCaloHit-earliestCaloHit;

      // Count all the tracker (Geiger) hits
      if (calData.has_calibrated_tracker_hits())
      {
        const snemo::datamodel::calibrated_data::tracker_hit_collection_type& trackerHits = calData.calibrated_tracker_hits();
        for (snemo::datamodel::calibrated_data::tracker_hit_collection_type::const_iterator   iHit = trackerHits.begin(); iHit != trackerHits.end(); ++iHit) {
          geigerHitCount++;
          // Get the location of the tracker hit
          const snemo::datamodel::calibrated_tracker_hit & hit = iHit->get();
          // Encode it into an integer so we can easily put it in an ntuple branch
          validation_.t_cell_hit_count_.push_back(EncodeLocation (hit));
        }
      }
    }
  catch (std::logic_error& e) {
    std::cerr << "failed to grab CD bank : " << e.what() << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  }

  // Number of tracker clusters comes from the TCD databank
  // We want two clusters of three cells
  try {
    const snemo::datamodel::tracker_clustering_data& clusterData = workItem.get<snemo::datamodel::tracker_clustering_data>("TCD");
    if (clusterData.has_default_solution ()) // Looks as if there is a possibility of alternative solutions. Is it sufficient to use the default?
      {
        snemo::datamodel::tracker_clustering_solution solution = clusterData.get_default_solution () ;
        snemo::datamodel::tracker_clustering_solution::cluster_col_type clusters=solution.get_clusters();
        for (snemo::datamodel::tracker_clustering_solution::cluster_col_type::const_iterator iCluster = clusters.begin();  iCluster != clusters.end(); ++ iCluster)
        {
          ++clusterCount;
        }
      }
  }
  catch (std::logic_error& e) {
    std::cerr << "failed to grab TCD bank : " << e.what() << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  }


  // Number of particle tracks PTD databank
  // We want two particle tracks to calculate 2e internal/external probability
  // If we have one track and a remote hit, we can calculate 1e1gamma probabilities

  try
  {
    const snemo::datamodel::particle_track_data& trackData = workItem.get<snemo::datamodel::particle_track_data>("PTD");
    if (trackData.has_particles ())
    {
      for (uint iParticle=0;iParticle<trackData.get_number_of_particles();++iParticle)
      {

        snemo::datamodel::particle_track track=trackData.get_particle(iParticle);
        switch (track.get_charge())
        {
            // Track count, positive and negative track count
          case snemo::datamodel::particle_track::NEGATIVE:
          {
            negativeTrackCount++;
            trackCount++;
            break;
          }
          case snemo::datamodel::particle_track::POSITIVE:
          {
            positiveTrackCount++;
            trackCount++;
            break;
          }
          case snemo::datamodel::particle_track::UNDEFINED:
          {
            trackCount++;
            break;
          }
          case snemo::datamodel::particle_track::NEUTRAL:
          {
            for (unsigned int hit=0; hit<track.get_associated_calorimeter_hits().size();++hit)
            {
              const snemo::datamodel::calibrated_calorimeter_hit & calo_hit = track.get_associated_calorimeter_hits().at(hit).get();
              double energy=calo_hit.get_energy();
              unassociatedEnergy+=energy;
              if (energy > LOW_ENERGY_LIMIT)
              {
                unassocOverThreshold += energy;
              }
            }
            continue;
          }
          default:
            continue;
        }
        // Number of tracker hits
        const snemo::datamodel::tracker_trajectory & the_trajectory = track.get_trajectory();
        const snemo::datamodel::tracker_cluster & the_cluster = the_trajectory.get_cluster();
        int numHits = the_cluster.get_number_of_hits();
        if (numHits>0)
        {
          allTrackHitCounts.push_back(numHits); // Vector of hits per track
          // Is the track associated to a calorimeter hit?
          if (track.has_associated_calorimeter_hits())
          {
            associatedTrackCount++;
          }
        }
      }
    }
  }// end try on PTD bank
  catch (std::logic_error& e) {
    std::cerr << "failed to grab PTD bank : " << e.what() << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  } //end catch

  validation_.h_total_calorimeter_energy_ = totalCalorimeterEnergy;
  validation_.h_calo_energy_over_threshold_ = energyOverThreshold;

  // Unassociated calorimeter energy is the total energy of the gammas

  validation_.h_unassociated_calorimeter_energy_ = unassociatedEnergy;
  validation_.h_unassociated_energy_over_threshold_ = unassocOverThreshold;
  validation_.h_associated_calorimeter_energy_ = totalCalorimeterEnergy-unassociatedEnergy;
  validation_.h_associated_energy_over_threshold_ = energyOverThreshold - unassocOverThreshold;

  // Timing
  validation_.h_calo_hit_time_separation_=TMath::Abs(timeDelay);

  // Counts
  validation_.h_calorimeter_hit_count_=caloHitCount;
  validation_.h_calo_hits_over_threshold_=nCalHitsOverLowLimit;
  validation_.h_geiger_hit_count_=geigerHitCount;
  validation_.h_cluster_count_=clusterCount;
  validation_.h_track_count_=trackCount;
  validation_.h_associated_track_count_=associatedTrackCount;
  validation_.h_negative_track_count_=negativeTrackCount;
  validation_.h_positive_track_count_=positiveTrackCount;
  validation_.v_all_track_hit_counts_=allTrackHitCounts;
  
  tree_->Fill();
  // MUST return a status, see ref dpp::processing_status_flags_type
  return dpp::base_module::PROCESS_OK;
}


void ValidationModule::ResetVars()
{
  validation_.v_all_track_hit_counts_.clear();
  validation_.t_cell_hit_count_.clear();
  validation_.c_calorimeter_hit_map_.clear();
  validation_.cm_average_calorimeter_energy_.clear();
  validation_.tm_average_drift_radius_.clear();
}

int ValidationModule::EncodeLocation(const snemo::datamodel::calibrated_tracker_hit & hit)
{
  int encodedLocation=hit.get_layer() + 100 * hit.get_row(); // There are fewer than 100 layers so this is OK
  if (hit.get_side()==0) encodedLocation = (-1 * encodedLocation) - 1;
  // Negative numbers are Italy side, positive are France side
  // layers are 0 to 8 with 0 being at the source foil and 8 by the main wall
  // rows go from 0 (mountain) to 112 (tunnel)
  return encodedLocation;
}


string ValidationModule::EncodeLocation(const snemo::datamodel::calibrated_calorimeter_hit & hit)
{ // Returns the geomID as a string ready for decoding
  std::stringstream buffer;
  buffer << hit.get_geom_id();
  return buffer.str();
}


//! [ValidationModule::reset]
void ValidationModule::reset() {
  hfile_->cd();
  tree_->Write();
  hfile_->Close(); //
  std::cout << "In reset: finished conversion, file closed " << std::endl;

  // clean up
  delete hfile_;
  filename_output_ = "Validation.root";
  this->_set_initialized(false);

}

