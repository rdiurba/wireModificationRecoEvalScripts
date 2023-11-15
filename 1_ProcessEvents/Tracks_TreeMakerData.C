// /pnfs/uboone/data/uboone/reconstructed/prod_v08_00_00_18/data_extbnb_mcc9.1_v08_00_00_18/run1_reco2_C1/00/00/69/58/PhysicsRun-2016_7_24_11_54_40-0006958-00040_20160724T214243_ext_bnb_20160724T234647_merged_20181107T121731_optfilter_20181224T075125_reco1_postwcct_postdl_20181224T082906_20190723T182808_reco2.root


/// Psuedo code:
/*

  Find all the T0 from acpttrigtagger
  
  then grab the tracks that have an association through acpttrigtagger
  
  using these tracks get the associated caloritmetry data product

  then got through all the calo's trajectory points, mark it's x, then find the associated hit and store it's attributes 
  

 */


// Standard things to include
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h> 
// These are the includes to use "Root" things 
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

// These are the larsoft includes that let you
// have access to data-products and the event 
// details
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
//#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/fwd.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

//I'll need, calo, tracks, hits, anab::T0
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/OpDetWaveform.h"

//This way you can be lazy
using namespace art;
using namespace std;

const std::vector<double> spikes{0.333333,
    0.50, 0.516667,0.533333,0.55,0.566667,0.583333,
    0.60, 0.616667,0.633333,0.65,0.666667,0.683333,
    0.70, 0.716667,0.733333,0.75,0.766667,0.783333,
    0.80, 0.816667,0.833333,0.85,0.866667,0.883333,
    0.90, 0.916667,0.933333,0.95,0.966667,0.983333,
    1.0, 1.333333, 1.5, 1.666667,
    2.0, 2.333333, 2.5, 2.666667,
    3.0, 3.333333, 3.5, 3.666667,
    4.0, 4.333333, 4.5, 4.666667,
    5.0, 5.333333, 5.5, 5.666667,
    6.0, 6.333333, 6.5, 6.666667,
    7.0, 7.333333, 7.5, 7.666667,
    8.0, 8.333333, 8.5, 8.666667,
    9.0, 9.333333, 9.5, 9.666667,
    10.0,10.333333,10.5,10.666667,
    11.0,11.333333,11.5,11.666667,
    12.0,12.333333,12.5,12.666667,
    13.0,13.333333,13.5,13.666667,
    14.0,14.333333,14.5,14.666667,
    15.0,15.333333,15.5,15.666667,
    16.0,16.333333,16.5,16.666667,
    17.0,17.333333,17.5,17.666667,
    18.0,18.333333,18.5,18.666667,
    19.0,19.333333,19.5,19.666667
    };

bool not_spiky(const std::vector<double>& spikes, const double hit_sigma){
  for(auto const& spike : spikes) {
    if(abs(spike-hit_sigma)<0.00001) {
      return false;
    }
  }
  return true;
}

double FoldAngle(double theta) {
  double th = std::abs(theta);
  if(th>0.5*util::pi()) th = util::pi() - th;
  return th;
}

void Tracks_TreeMakerData(){

  // create a vector of files we want to process
  //std::vector<std::string> filenames;

  // read in a file list that we get from a sam-def but remember it 
  // is very long so if we want to run over it all it'll take a while
  // and we'll probably want to break it up on the grid
  
  //ifstream myfile("file_list_4b.txt");
  //copy(istream_iterator<string>(myfile),
  //     istream_iterator<string>(),
  //     back_inserter(filenames));
  

  //We'll just check the first 10k files for now (that's probably ~25k tracks)
  //filenames.erase(filenames.begin()+200,filenames.end());
  
  //Here I just hard coding a file for testing, we can adjust this later 
  std::vector<std::string> filenames {"test.root"};
  
  // Here we will create all of our histograms 
  // I did this crazy inefficiently but I don't really care
  // This is currently only set up for single dimensional 
  // projections but extenting it to 2D will be straight forward
  TFile* out = new TFile("output.root","RECREATE");  
  TTree* fTree = new TTree("hit_v_track","HitPropertiesTrackProperties");

  int Run;
  fTree->Branch("Run",&Run);
  int Subrun;
  fTree->Branch("Subrun",&Subrun);
  int Event;
  fTree->Branch("Event",&Event);
  int Ntrks;
  fTree->Branch("Ntrks",&Ntrks);
  int Nacpt;
  fTree->Branch("Nacpt",&Nacpt);  
  int n_acpt;
  fTree->Branch("n_acpt",&n_acpt);  
  int Nhits;
  fTree->Branch("Nhits",&Nhits);
  int Nhits_Plane0;
  int Nhits_Plane1;
  int Nhits_Plane2;
  fTree->Branch("Nhits_Plane0",&Nhits_Plane0);
  fTree->Branch("Nhits_Plane1",&Nhits_Plane1);
  fTree->Branch("Nhits_Plane2",&Nhits_Plane2);
  int n_hits;
  fTree->Branch("n_hits",&n_hits);  
  double hit_Q;
  double hit_A;
  double hit_sigma;
  double hit_time;
  int hit_plane;
  int hit_isNotSpiky;
  int hit_isMC;
  fTree->Branch("hit_Q",&hit_Q);  
  fTree->Branch("hit_A",&hit_A);  
  fTree->Branch("hit_sigma",&hit_sigma);  
  fTree->Branch("hit_time",&hit_time);  
  fTree->Branch("hit_plane",&hit_plane);  
  fTree->Branch("hit_isNotSpiky",&hit_isNotSpiky);  
  fTree->Branch("hit_isMC",&hit_isMC);  

  double trk_x;
  double trk_y;
  double trk_z;
  double trk_L;
  double trk_StartX;
  double trk_StartY;
  double trk_StartZ;
  double trk_EndX;
  double trk_EndY;
  double trk_EndZ;
  double trk_cosTheta;
  double trk_phi;
  double trk_ThetaXZ_Plane0;
  double trk_ThetaYZ_Plane0;
  double trk_ThetaXZ_Plane1;
  double trk_ThetaYZ_Plane1;
  double trk_ThetaXZ_Plane2;
  double trk_ThetaYZ_Plane2;
  double trk_dEdx;
  double trk_RR;
  double trk_fracMC;

  fTree->Branch("trk_x",&trk_x);  
  fTree->Branch("trk_y",&trk_y);  
  fTree->Branch("trk_z",&trk_z);  
  fTree->Branch("trk_L",&trk_L);  

  fTree->Branch("trk_StartX",&trk_StartX);  
  fTree->Branch("trk_StartY",&trk_StartY);  
  fTree->Branch("trk_StartZ",&trk_StartZ);  

  fTree->Branch("trk_EndX",&trk_EndX);  
  fTree->Branch("trk_EndY",&trk_EndY);  
  fTree->Branch("trk_EndZ",&trk_EndZ);  

  fTree->Branch("trk_cosTheta",&trk_cosTheta);  
  fTree->Branch("trk_phi",&trk_phi);  
  fTree->Branch("trk_ThetaXZ_Plane0",&trk_ThetaXZ_Plane0);  
  fTree->Branch("trk_ThetaYZ_Plane0",&trk_ThetaYZ_Plane0);  
  fTree->Branch("trk_ThetaXZ_Plane1",&trk_ThetaXZ_Plane1);  
  fTree->Branch("trk_ThetaYZ_Plane1",&trk_ThetaYZ_Plane1);  
  fTree->Branch("trk_ThetaXZ_Plane2",&trk_ThetaXZ_Plane2);  
  fTree->Branch("trk_ThetaYZ_Plane2",&trk_ThetaYZ_Plane2);  
  fTree->Branch("trk_dEdx",&trk_dEdx);  
  fTree->Branch("trk_RR",&trk_RR);  
  fTree->Branch("trk_fracMC",&trk_fracMC);  


  //Addition from Wes: Wire metrics
  int channel;
  double wire_Q;
  double wire_A;
  double wire_sigma;
  double wire_time;
  int wire_begin;
  int wire_end;
  double subwire_Q;
  double subwire_A;
  double subwire_sigma;
  double subwire_time;
  int subwire_begin;
  int subwire_end;

  fTree->Branch("channel",&channel);
  fTree->Branch("wire_Q",&wire_Q);
  fTree->Branch("wire_A",&wire_A);
  fTree->Branch("wire_sigma",&wire_sigma);
  fTree->Branch("wire_time",&wire_time);
  fTree->Branch("wire_begin",&wire_begin);
  fTree->Branch("wire_end",&wire_end);
  fTree->Branch("subwire_Q",&subwire_Q);
  fTree->Branch("subwire_A",&subwire_A);
  fTree->Branch("subwire_sigma",&subwire_sigma);
  fTree->Branch("subwire_time",&subwire_time);
  fTree->Branch("subwire_begin",&subwire_begin);
  fTree->Branch("subwire_end",&subwire_end);
  
  //First things fist we need to iterate through each event
  // gallery makes it easy to just hand a vector of files
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    
    /// Prep our Branches
    Run = 0; //
    Subrun = 0; //
    Event = 0; //
    Ntrks = 0; //
    Nacpt = 0;//
    n_acpt = -1;//
    Nhits = 0;//
    Nhits_Plane0 = 0;//
    Nhits_Plane1 = 0;//
    Nhits_Plane2 = 0;//
    n_hits = -1;//
    hit_Q = 0; //
    hit_A = 0; //
    hit_sigma = 0; //
    hit_time = 0; //
    hit_plane = 0; //
    hit_isNotSpiky = -1; //
    trk_x = 0;//
    trk_y = 0;//
    trk_z = 0;//
    trk_cosTheta = 0;//
    trk_phi = 0;//
    trk_ThetaXZ_Plane0 = 0;//
    trk_ThetaYZ_Plane0 = 0;//
    trk_ThetaXZ_Plane1 = 0;//
    trk_ThetaYZ_Plane1 = 0;//
    trk_ThetaXZ_Plane2 = 0;//
    trk_ThetaYZ_Plane2 = 0;//
    trk_dEdx = 0;//
    trk_RR = 0;

    trk_StartX = 0;//
    trk_StartY = 0;//
    trk_StartZ = 0;//
    trk_EndX = 0;//
    trk_EndY = 0;//
    trk_EndZ = 0; //
    
    
    trk_L = 0;//
    trk_fracMC = 0;//
    hit_isMC = 0;//

    //wire stuff
    channel = -1;
    wire_Q = 0;
    wire_A = 0;
    wire_sigma = 0;
    wire_time = 0;
    wire_begin = -1;
    wire_end = -1;    
    subwire_Q = 0;
    subwire_A = 0;
    subwire_sigma = 0;
    subwire_time = 0;
    subwire_begin = -1;
    subwire_end = -1;    
    

    Run = ev.eventAuxiliary().run();
    Subrun = ev.eventAuxiliary().subRun();
    Event = ev.eventAuxiliary().event();
    
       
    //We will now found the events that have a ACPT in-time track
    auto const& t0s = *ev.getValidHandle<vector<anab::T0>>("acpttrigtagger");
    
    //Skipping those that don't
    if(t0s.size() == 0) continue;
    
    Nacpt = t0s.size();
    
    //This associates the T0 tag with the track
    auto const &t0_assoc_handle =
      ev.getValidHandle<art::Assns<anab::T0, recob::Track>>("acpttrigtagger");

    //Make a vector that will hold our tagged tracks
    std::vector<recob::Track> ACPT_tracks;
    
    // store the tagged tracks into that vector
    for(auto &ass : *t0_assoc_handle){
      art::Ptr<recob::Track> temp_trk = ass.second;
      ACPT_tracks.emplace_back(*temp_trk);
    }
    
    // Now we'll need to set things up to collect the calorimetry data
    // Start by saying which tracks we want:
    auto const & track_list_ptr = ev.getValidHandle<std::vector <recob::Track> >("pandora");
    auto const & track_list = (*track_list_ptr);
    
    Ntrks = track_list.size();

    //This let's us find which hits are assocaited to a give track trajectory point
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(track_list_ptr, ev, "pandora"); 

    //get hit assn to recob::Wire
    auto const& hit_handle = ev.getValidHandle<std::vector<recob::Hit>>("gaushit");
    art::FindMany<recob::Wire> wires_per_hit(hit_handle,ev,"gaushit");
    

    //Let's loop through our tracks and find our calorimetry things
    for(auto &trk : ACPT_tracks){      
      n_acpt++;
      for(int itrk = 0; itrk < int(track_list.size()); itrk++){	
	if(trk.ID() == track_list.at(itrk).ID()){ 
	  
	  // Now we have a track which is matched to a ACPT crossing track
	  // now we want the hits for this track

	  // This is the vector of hits
	  auto vhit = fmthm.at(itrk);
	  Nhits = vhit.size();

	  // Calculate number of non-spiky hits per plane
	  Nhits_Plane0 = 0;
	  Nhits_Plane1 = 0;
	  Nhits_Plane2 = 0;
	  for( unsigned int i_h = 0; i_h < vhit.size(); i_h++ ) {
	    if( not( not_spiky(spikes, vhit[i_h]->RMS()) ) ) continue;
	    if( vhit[i_h]->WireID().Plane == 0 ) Nhits_Plane0++;
	    if( vhit[i_h]->WireID().Plane == 1 ) Nhits_Plane1++;
	    if( vhit[i_h]->WireID().Plane == 2 ) Nhits_Plane2++;
	  }

	  //Next we'll get the associated calorimetries things:
	  art::FindMany<anab::Calorimetry>  fmcal(track_list_ptr, ev, "pandoracaliSCE");//"pandoracali");
	  
 	  trk_L = trk.Length();
	  trk_StartX = trk.Start().X();
	  trk_StartY = trk.Start().Y();
	  trk_StartZ = trk.Start().Z();
	  trk_EndX = trk.End().X();
	  trk_EndY = trk.End().Y();
	  trk_EndZ = trk.End().Z();


	  // This is the vector of traj point info
	  // the Index() of this is the traj point of this track
	  auto vmeta = fmthm.data(itrk);
	  
	  // Now we can get the calorimetry points
	  std::vector<const anab::Calorimetry*> calos = fmcal.at(itrk);

	  // this will count which calo point we're on
	  int count = 0;

	  //iterate through the planes:
	  for(int pl = 0; pl < 3; pl++){
	    	 
	    //iterate through track meta points :	    
	    for(int vp = 0; vp < vmeta.size(); vp++){
	      
	      // store the track trajectory point index
	      // for the track meta point
	      int ind = vmeta[vp]->Index();

	      // check that the traj point is in the calorimetry point
	      // and belongs to the plane we are interested in 
	      if(track_list.at(itrk).HasValidPoint(ind) && vhit[vp]->WireID().Plane == pl){
		
		n_hits++;

		// Grab the track traj point
		// WE DON'T CURRENTLY USE THIS
		// I kept it for testing purposes 
		auto trjp = track_list.at(itrk).TrajectoryPoint(ind);
		
		// Grab the calo point
		auto calp = calos[pl]->XYZ()[count];
		auto caldEdx = calos[pl]->dEdx()[count];
		auto calRR = calos[pl]->ResidualRange()[count];
		
		// We need to calculate the angles 
		// of the calo points 
		double Phi = 0;
		double cosTheta = 0;
		double ThetaXZ_Plane0 = 0;
		double ThetaYZ_Plane0 = 0;
		double ThetaXZ_Plane1 = 0;
		double ThetaYZ_Plane1 = 0;
		double ThetaXZ_Plane2 = 0;
		double ThetaYZ_Plane2 = 0;
		
		if(count < vmeta.size()-1){
		  auto angle = (calos[pl]->XYZ()[count]) - (calos[pl]->XYZ()[count+1]);
		  Phi = angle.Phi(); 
		  cosTheta = cos(angle.Theta());
		  ThetaXZ_Plane2 = atan2(angle.X(),angle.Z());
		  ThetaYZ_Plane2 = atan2(angle.Y(),angle.Z());		  
		  // Calculate detector angles on other planes
		  ThetaXZ_Plane0 = atan2( angle.X(), -std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaYZ_Plane0 = atan2( 0.5*angle.Y()+std::sqrt(0.75)*angle.Z(), -std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaXZ_Plane1 = atan2( angle.X(), std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaYZ_Plane1 = atan2( 0.5*angle.Y()-std::sqrt(0.75)*angle.Z(), std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		}
		else{
		  auto angle = (calos[pl]->XYZ()[count-1]) - (calos[pl]->XYZ()[count]);  
		  Phi = angle.Phi(); 
		  cosTheta = cos(angle.Theta());
		  ThetaXZ_Plane2 = atan2(angle.X(),angle.Z());
		  ThetaYZ_Plane2 = atan2(angle.Y(),angle.Z());		  
		  // Calculate detector angles on other planes
		  ThetaXZ_Plane0 = atan2( angle.X(), -std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaYZ_Plane0 = atan2( 0.5*angle.Y()+std::sqrt(0.75)*angle.Z(), -std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaXZ_Plane1 = atan2( angle.X(), std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		  ThetaYZ_Plane1 = atan2( 0.5*angle.Y()-std::sqrt(0.75)*angle.Z(), std::sqrt(0.75)*angle.Y()+0.5*angle.Z() );
		}
		// Fold all angles
		ThetaXZ_Plane0 = FoldAngle(ThetaXZ_Plane0);
		ThetaYZ_Plane0 = FoldAngle(ThetaYZ_Plane0);
		ThetaXZ_Plane1 = FoldAngle(ThetaXZ_Plane1);
		ThetaYZ_Plane1 = FoldAngle(ThetaYZ_Plane1);
		ThetaXZ_Plane2 = FoldAngle(ThetaXZ_Plane2);
		ThetaYZ_Plane2 = FoldAngle(ThetaYZ_Plane2);

		// Grab the matched hit
		auto hit = vhit[vp];


		hit_A = hit->PeakAmplitude();
		hit_Q = hit->Integral();
		hit_sigma = hit->RMS();
		hit_time = hit->PeakTime();
		hit_plane = pl;
		hit_isNotSpiky = (int)not_spiky(spikes, hit_sigma);
		
		trk_x = calp.X(); 
		trk_y = calp.Y(); 
		trk_z = calp.Z(); 
		trk_phi = Phi;
		trk_cosTheta = cosTheta;
		trk_dEdx = caldEdx;
		trk_RR   = calRR;
		trk_ThetaXZ_Plane0 = ThetaXZ_Plane0;
		trk_ThetaYZ_Plane0 = ThetaYZ_Plane0;
		trk_ThetaXZ_Plane1 = ThetaXZ_Plane1;
		trk_ThetaYZ_Plane1 = ThetaYZ_Plane1;
		trk_ThetaXZ_Plane2 = ThetaXZ_Plane2;
		trk_ThetaYZ_Plane2 = ThetaYZ_Plane2;

		//wire info...
		channel = hit->Channel();
                
		//Grab wire
		auto wire_vec = wires_per_hit.at(hit.key());
		
		//there should be only one associate recob::Wire...
		if(wire_vec.size()==0){
		  std::cout << "ERROR: No associated wire!" << std::endl;
		  continue;
		}
		else if(wire_vec.size()>1){
		  std::cout << "WARNING: More than one associated Wire? Only taking the first." << std::endl;
		}
		
		auto const* wire_ptr = wire_vec[0];
		auto const& signal_rois = wire_ptr->SignalROI();
		for (auto const& range: signal_rois.get_ranges()) {
		  if(hit->PeakTime()<range.begin_index() || 
		     hit->PeakTime()>range.begin_index()+range.size()) continue;
		  wire_begin = range.begin_index();
		  wire_end = wire_begin+range.size();
		  wire_Q = 0;
		  wire_A = 0;
		  wire_time = 0;
		  wire_sigma = 0;
		  
		  subwire_begin = hit->PeakTime()-hit->RMS();
		  if(subwire_begin<range.begin_index()) subwire_begin = range.begin_index();
		  subwire_end = hit->PeakTime()+hit->RMS();
		  if(subwire_end==subwire_begin) subwire_end+=1;
		  if(subwire_end>=range.begin_index()+range.size()) subwire_end = range.begin_index()+range.size();
		  
		  subwire_Q = 0;
		  subwire_A = 0;
		  subwire_time = 0;
		  subwire_sigma = 0;
		  
		  for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
		    auto val = range.data().at(i_t);
		    auto tick = wire_begin+i_t;
		    if(val>wire_A) wire_A = val;
		    wire_Q += val;
		    wire_time += val*tick;

		    if( tick>=subwire_begin && tick<subwire_end ){
		      if(val>subwire_A) subwire_A = val;
		      subwire_Q += val;
		      subwire_time += val*tick;
		    }

		  }
		  wire_time = wire_time/wire_Q;
		  subwire_time = subwire_time/subwire_Q;
		  for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
		    auto val = range.data().at(i_t);
		    auto tick = wire_begin+i_t;
		    wire_sigma += val*(tick-wire_time)*(tick-wire_time);

		    if( tick>=subwire_begin && tick<subwire_end )
		      subwire_sigma += val*(tick-subwire_time)*(tick-subwire_time);
		  }
		  wire_sigma = std::sqrt(wire_sigma/wire_Q);
		  subwire_sigma = std::sqrt(subwire_sigma/subwire_Q);
		}



		fTree->Fill();
		
		// this tracks the correct calorimetry 
		// we are supposed to be anlayzing 
		count++;
	      }
	    }
	    count = 0; 
	  } // end loop over planes
	}
      }
    } // end loop over ACPT tracks
  } // end loop over files

  out->cd();  
  fTree->Write();
  

}
  

