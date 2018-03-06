#include "ValidationModule.h"
const double highEnergyLimit=0.150;// 150 keV
const double lowEnergyLimit=0.050; // 50 keV
const double electronMass=0.5109989461; // From pdg
const double speedOfLight=299792458 * 1e-9 * 1000; // Millimeters per nanosecond

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
  tree_->Branch("reco.calorimeter_hit_count",&Validation_.calorimeter_hit_count_);
  tree_->Branch("reco.cluster_count",&Validation_.cluster_count_);
  tree_->Branch("reco.track_count",&Validation_.track_count_);
  tree_->Branch("reco.negative_track_count",&Validation_.negative_track_count_);
  tree_->Branch("reco.positive_track_count",&Validation_.positive_track_count_);
  tree_->Branch("reco.associated_track_count",&Validation_.associated_track_count_);
  tree_->Branch("reco.geiger_hit_count",&Validation_.geiger_hit_count_);
  tree_->Branch("reco.all_track_hit_counts",&Validation_.all_track_hit_counts_);
  
  // Energies and calo times
  tree_->Branch("reco.total_calorimeter_energy",&Validation_.total_calorimeter_energy_);
  tree_->Branch("reco.unassociated_calorimeter_energy",&Validation_.unassociated_calorimeter_energy_);
  tree_->Branch("reco.calo_hit_time_separation",&Validation_.calo_hit_time_separation_);
  
  this->_set_initialized(true);
}
//! [ValidationModule::Process]
dpp::base_module::process_status
ValidationModule::process(datatools::things& workItem) {
  
  
  // declare internal variables to mimic the ntuple variables, names are same but in camel case
  double totalCalorimeterEnergy=0;
  double timeDelay=-1;
  int clusterCount=0;
  int trackCount=0;
  int alphaCount=0;
  int associatedTrackCount=0;
  int caloHitCount=0;
  double earliestCaloHit=-1.;
  double latestCaloHit=-1.;
  int geigerHitCount=0;
  int negativeTrackCount=0;
  int positiveTrackCount=0;
  std::vector<int> allTrackHitCounts;

  // Grab calibrated data bank
  // Calibrated data will only be present in reconstructed files,
  // so wrap in a try block
  
  try {
    const snemo::datamodel::calibrated_data& calData = workItem.get<snemo::datamodel::calibrated_data>("CD");

    int nCalorimeterHits=0;
    int nCalHitsOverHighLimit=0;
    int nCalHitsOverLowLimit=0;
    
    if (calData.has_calibrated_calorimeter_hits())
      {
        const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & calHits=calData.calibrated_calorimeter_hits();
        for (snemo::datamodel::calibrated_data::calorimeter_hit_collection_type::const_iterator   iHit = calHits.begin(); iHit != calHits.end(); ++iHit) {
          const snemo::datamodel::calibrated_calorimeter_hit & calHit = iHit->get();
          double energy=calHit.get_energy() ;
          totalCalorimeterEnergy += energy;

          double hitTime=calHit.get_time();
          if (hitTime > latestCaloHit)
          {
            latestCaloHit = hitTime;
          }
          if (hitTime < earliestCaloHit || earliestCaloHit <= 0 )
          {
            earliestCaloHit=hitTime;
          }
          
          ++nCalorimeterHits;
          if (energy>=highEnergyLimit)++nCalHitsOverHighLimit;
          if (energy>=lowEnergyLimit)++nCalHitsOverLowLimit;
        }
      }
    
    timeDelay=latestCaloHit-earliestCaloHit;
    caloHitCount=nCalHitsOverLowLimit;

      // Count all the tracker (Geiger) hits
      if (calData.has_calibrated_tracker_hits())
      {
        const snemo::datamodel::calibrated_data::tracker_hit_collection_type& trackerHits = calData.calibrated_tracker_hits();
        for (snemo::datamodel::calibrated_data::tracker_hit_collection_type::const_iterator   iHit = trackerHits.begin(); iHit != trackerHits.end(); ++iHit) {
          geigerHitCount++;
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
        
        // Basic debug info for any track
        if (track.get_charge()==snemo::datamodel::particle_track::NEGATIVE)
        {
          negativeTrackCount++;
          trackCount++;
        }
        if (track.get_charge()==snemo::datamodel::particle_track::POSITIVE)
        {
          positiveTrackCount++;
          trackCount++;
        }
        if (track.get_charge()==snemo::datamodel::particle_track::UNDEFINED)
        {
          trackCount++;
        }
        else continue;
      }
    }
  }// end try on PTD bank
  catch (std::logic_error& e) {
    std::cerr << "failed to grab PTD bank : " << e.what() << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  } //end catch

  
  // Initialise variables that might not otherwise get set
  // It does not restart the vector for each entry so we have to do that manually
  ResetVars();

  Validation_.total_calorimeter_energy_ = totalCalorimeterEnergy;

  // Unassociated calorimeter energy is the total energy of the gammas
  double unassociatedEnergy=0.;
//  for (int i=0; i<gammaEnergies.size();i++)
//  {
//    unassociatedEnergy += gammaEnergies.at(i);
//  }
//
  Validation_.unassociated_calorimeter_energy_ = unassociatedEnergy;

  // Timing
  Validation_.calo_hit_time_separation_=TMath::Abs(timeDelay);

  // Debug information
  Validation_.calorimeter_hit_count_=caloHitCount;
  Validation_.geiger_hit_count_=geigerHitCount;
  Validation_.cluster_count_=clusterCount;
  Validation_.track_count_=trackCount;
  Validation_.negative_track_count_=negativeTrackCount;
  Validation_.positive_track_count_=positiveTrackCount;
  //Validation_.associated_track_count_=electronCandidates.size();
  
  tree_->Fill();
  
  // MUST return a status, see ref dpp::processing_status_flags_type
  return dpp::base_module::PROCESS_OK;
}


void ValidationModule::ResetVars()
{

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

