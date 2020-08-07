// -*- C++ -*-
//
// Package:    MyAna/L1TPAnalyzer
// Class:      L1TPAnalyzer
//
/**\class L1TPAnalyzer L1TPAnalyzer.cc MyAna/L1TPAnalyzer/plugins/L1TPAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Brandon Chiarito
//         Created:  Wed, 13 May 2020 23:47:18 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

// ROOT
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibCalorimetry/CaloTPG/src/CaloTPGTranscoderULUT.h"

//#include <algorithm>
//#include "TVector3.h"
//#include "TLorentzVector.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TTree.h"
//#include "TMath.h"

//samplingFactor for ietaAbs 1 to 16. Source: https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/SimCalorimetry/HcalSimProducers/python/hcalSimParameters_cfi.py#L60
std::vector <float> samplingFactors_hb = {125.44, 125.54, 125.32, 125.13, 124.46, 125.01, 125.22, 125.48, 124.45, 125.90, 125.83, 127.01, 126.82, 129.73, 131.83, 143.52};
//samplingFactor for ietaAbs 16 to 29. Source: https://github.com/cms-sw/cmssw/blob/CMSSW_10_2_X/SimCalorimetry/HcalSimProducers/python/hcalSimParameters_cfi.py#L77
std::vector <float> samplingFactors_he = {210.55, 197.93, 186.12, 189.64, 189.63, 190.28, 189.61, 189.60, 190.12, 191.22, 190.90, 193.06, 188.42, 188.42};

using std::vector;

class L1TPAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1TPAnalyzer(const edm::ParameterSet&);
      ~L1TPAnalyzer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      bool cPrintAllEvents;
      bool cPrintCSVFile;
      int cTowerIeta;
      int cTowerIphi;

      int my_event_label;
      edm::InputTag emulTPsTag_;
      edm::EDGetTokenT<HcalTrigPrimDigiCollection> tok_emulTPs_;
      edm::EDGetTokenT<std::vector<PCaloHit>> simhits_;
  
      // Main Ntuple Ttree and branches
      TTree *fTree;
      Int_t fmyev;
      Int_t fEventNum;
      Int_t fRunNum;
      Int_t fLumiNum;
      vector<Double_t> fiphi;
      vector<Double_t> fieta;
      vector<Double_t> feta;
      vector<Double_t> fphi;
      vector<Double_t> feta_t;
      vector<Double_t> fphi_t;
      vector<Double_t> fdepth;
      vector<Double_t> fsubdet;
      vector<Double_t> fversion;
      vector<Double_t> frawid;
      vector<Double_t> fsoi;
      vector<Double_t> fsize;
      vector<Double_t> ffg_all;
      vector<Double_t> fEt_soi;
      vector<Double_t> fEt_soi_rev;
      vector<Double_t> fEt_s0;
      vector<Double_t> fEt_s1;
      vector<Double_t> fEt_s2;
      vector<Double_t> fEt_s3;

      vector<Double_t> fsimhits_rawid;
      vector<Double_t> fsimhits_time;
      vector<Double_t> fsimhits_energy;
      vector<Double_t> fsimhits_energy_coded;
      vector<Double_t> fsimhits_EMfrac;
      vector<Double_t> fsimhits_bx;
      vector<Double_t> fsimhits_event;
      vector<Double_t> fsimhits_ieta;
      vector<Double_t> fsimhits_iphi;
      vector<Double_t> fsimhits_eta;
      vector<Double_t> fsimhits_phi;
      vector<Double_t> fsimhits_depth;
      vector<Double_t> fsimhits_subdet;

      TH2D * ftp_2d;
      TH2D * ftp_2d_linear;
      TH2D * fsimhit_2d;
      TH2D * fsimhit_2d_coded;
      TGraph * flinear_compare;
      TGraph * fcompressed_compare;
      TH2D * flinear_compare_hist;
      TH2D * fcompressed_compare_hist;
};

L1TPAnalyzer::L1TPAnalyzer(const edm::ParameterSet& iConfig)
 :
  cPrintAllEvents(iConfig.getUntrackedParameter<bool>("printAllEvents")),
  cPrintCSVFile(iConfig.getUntrackedParameter<bool>("printCSVFile")),
  cTowerIeta(iConfig.getUntrackedParameter<int>("towerIeta")),
  cTowerIphi(iConfig.getUntrackedParameter<int>("towerIphi"))
{
    emulTPsTag_ = iConfig.getParameter<edm::InputTag>("emulTPs");
    my_event_label = 0;

    tok_emulTPs_ = consumes<HcalTrigPrimDigiCollection>(emulTPsTag_);
    simhits_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits","HcalHits","SIM"));

    edm::Service<TFileService> fs;
    fTree = fs->make<TTree>("fTree","fTree");
    fTree->Branch("myev",&fmyev,"myev/I");
    fTree->Branch("eventNum",&fEventNum,"eventNum/I");
    fTree->Branch("runNum",&fRunNum,"runNum/I");
    fTree->Branch("lumiNum",&fLumiNum,"lumiNum/I");
    fTree->Branch("iphi",&fiphi);
    fTree->Branch("ieta",&fieta);
    fTree->Branch("eta",&feta);
    fTree->Branch("phi",&fphi);
    fTree->Branch("eta_t",&feta_t);
    fTree->Branch("phi_t",&fphi_t);
    fTree->Branch("depth",&fdepth);
    fTree->Branch("subdet",&fsubdet);
    fTree->Branch("version",&fversion);
    fTree->Branch("rawid",&frawid);
    fTree->Branch("soi",&fsoi);
    fTree->Branch("size",&fsize);
    fTree->Branch("fg_all",&ffg_all);
    fTree->Branch("Et_soi",&fEt_soi);
    fTree->Branch("Et_soi_rev",&fEt_soi_rev);
    fTree->Branch("Et_s0",&fEt_s0);
    fTree->Branch("Et_s1",&fEt_s1);
    fTree->Branch("Et_s2",&fEt_s2);
    fTree->Branch("Et_s3",&fEt_s3);

    fTree->Branch("simhits_rawid",&fsimhits_rawid);
    fTree->Branch("simhits_time",&fsimhits_time);
    fTree->Branch("simhits_energy",&fsimhits_energy);
    fTree->Branch("simhits_energy_coded",&fsimhits_energy_coded);
    fTree->Branch("simhits_EMfrac",&fsimhits_EMfrac);
    fTree->Branch("simhits_bx",&fsimhits_bx);
    fTree->Branch("simhits_event",&fsimhits_bx);
    fTree->Branch("simhits_ieta",&fsimhits_ieta);
    fTree->Branch("simhits_iphi",&fsimhits_iphi);
    fTree->Branch("simhits_eta",&fsimhits_eta);
    fTree->Branch("simhits_phi",&fsimhits_phi);
    fTree->Branch("simhits_depth",&fsimhits_depth);
    fTree->Branch("simhits_subdet",&fsimhits_subdet);

    ftp_2d = fs->make<TH2D>("tp_2d","tp_2d",100,-50,50,100,0,100);
    ftp_2d_linear = fs->make<TH2D>("tp_2d_linear","tp_2d_linear",100,-50,50,100,0,100);
    fsimhit_2d = fs->make<TH2D>("simhit_2d","simhit_2d",100,-50,50,100,0,100);
    fsimhit_2d_coded = fs->make<TH2D>("simhit_2d_coded","simhit_2d_coded",100,-50,50,100,0,100);

    flinear_compare = fs->make<TGraph>(1);
    flinear_compare->SetName("linear_compare");
    fcompressed_compare = fs->make<TGraph>(1);
    fcompressed_compare->SetName("compressed_compare");

    int compare_binsx = 300;
    int compare_xlow = 0;
    int compare_xhigh = 1000;
    int compare_binsy = 300;
    int compare_ylow = 0;
    int compare_yhigh = 1000;
    flinear_compare_hist = fs->make<TH2D>("linear_compare_hist","linear_compare_hist",compare_binsx,compare_xlow,compare_xhigh,compare_binsy,compare_ylow,compare_yhigh);
    fcompressed_compare_hist = fs->make<TH2D>("compressed_compare_hist","compressed_compare_hist",compare_binsx,compare_xlow,compare_xhigh,compare_binsy,compare_ylow,compare_yhigh);
}


L1TPAnalyzer::~L1TPAnalyzer()
{
}




// ------------ method called for each event  ------------
void
L1TPAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    my_event_label += 1;
    fmyev = my_event_label;

    using namespace edm;
    using std::cout;
    using std::endl;
    using std::vector;

    if(!cPrintCSVFile && cPrintAllEvents) cout << "\n\nevent " << iEvent.id().event() << " lumi " <<  iEvent.id().luminosityBlock() << " run " << iEvent.id().run() << "\n" << endl;

    fEventNum = iEvent.id().event();
    fRunNum = iEvent.id().run();
    fLumiNum = iEvent.id().luminosityBlock();

    edm::Handle<HcalTrigPrimDigiCollection> emulTPs;
    iEvent.getByToken(tok_emulTPs_, emulTPs);

    edm::Handle<std::vector<PCaloHit>> simhits;
    iEvent.getByToken(simhits_, simhits);

    ESHandle<HcalDDDRecConstants> pHRNDC;
    iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
    const HcalDDDRecConstants *hcons = &(*pHRNDC);

    edm::ESHandle<HcalTrigTowerGeometry> trigGeometry;
    iSetup.get<CaloGeometryRecord>().get(trigGeometry);
    const HcalTrigTowerGeometry *tgeo = &(*trigGeometry);

    // outcoder
    edm::ESHandle<CaloTPGTranscoder> outTranscoder;
    iSetup.get<CaloTPGRecord>().get(outTranscoder);
    const HcalTPGCompressor * outcoder_ = outTranscoder->getHcalCompressor().get();

    const CaloTPGTranscoderULUT * outTranscoder_cast = (const CaloTPGTranscoderULUT *) &(*outTranscoder);

    // build outcoder map
    
//        int compressed_energy = (outcoder_->compress(ids[0], int(simhit.energy() * samplingFactor), false)).compressedEt();
/*
    int range1 = 100;
    int range2 = 100;
    int range3 = 200;
    vector<vector<vector<vector<int>>>> outcoder_revmap( range1, vector<vector<vector<int>>>( range2, vector<vector<int>>( range3, vector<int>() ) ) );

    for(int ie = 1; ie < 100; ie++) {
      for(int ip = 1; ip < 80; ip++) {
        for(int r = 0; r < 1000; r++) {
          cout << "ie " << ie << " ip " << ip << " r " << r;
          int outcoder_from = r;
          HcalTrigTowerDetId trigtower_ieip(ie, ip);
          int outcoder_to = int( (outcoder_->compress(trigtower_ieip, outcoder_from, false)).compressedEt() );
          cout << " result " << outcoder_to << endl;
          outcoder_revmap[ie][ip][outcoder_to].push_back(outcoder_from);
        }
      }
    }*/
    //cout << outcoder_revmap.at(1).at(1).at(1) << endl;


   /* std::map<std::tuple<int,int,int>, vector<int>> outcoder_revmap;

    for(int ie = -28; ie < 28; ie++) {
      if(ie == 0) continue;
      cout << "ie " << ie << endl;
      for(int ip = 1; ip < 80; ip++) {
        if(ip == 0) continue;
        cout << "  ip " << ip << endl;
        for(int r = 0; r < 2; r++) {
          cout << "    r " << r << endl;

          //cout << "ie " << ie << " ip " << ip << " r " << r;
          int outcoder_from = r;
          HcalTrigTowerDetId trigtower_ieip(ie, ip);
          //HcalTrigTowerDetId * trigtower_ieip = new HcalTrigTowerDetId(ie, ip);
          //delete trigtower_ieip;

          cout << "    going to call compress()" << endl;
          int outcoder_to = int( (outcoder_->compress(trigtower_ieip, outcoder_from, false)).compressedEt() );
          outcoder_to = 10;
          cout << "      result " << outcoder_to << endl;

          //int outcoder_to = int( (outcoder_->compress(HcalTrigTowerDetId(ie,ip), outcoder_from, false)).compressedEt() );

          //outcoder_->compress(*trigtower_ieip, outcoder_from, false);
          //int outcoder_to = int( result.compressedEt() );
          

          //struct argument arg;
          //arg.ieta = ie;
          //arg.iphi = ip;
          //arg.to_value = outcoder_to;

          //std::tuple<int,int,int> arg = std::make_tuple(ie, ip, outcoder_to);

          //if(outcoder_revmap.find( arg ) == outcoder_revmap.end()) {
            //cout << "  new vector" << endl;
            //outcoder_revmap[arg] = vector<int>();
            //outcoder_revmap[arg].push_back(outcoder_to);
          //} else {
            //cout << "  add to existting vector" << endl;
            //outcoder_revmap[arg].push_back(outcoder_to);
          //}
          cout << "      done" << endl;

        }
      }
    }*/

    fiphi.clear();
    fieta.clear();
    feta.clear();
    fphi.clear();
    feta_t.clear();
    fphi_t.clear();
    fdepth.clear();
    fsubdet.clear();
    fversion.clear();
    frawid.clear();
    fsoi.clear();
    fsize.clear();
    ffg_all.clear();
    fEt_soi.clear();
    fEt_soi_rev.clear();
    fEt_s0.clear();
    fEt_s1.clear();
    fEt_s2.clear();
    fEt_s3.clear();

    fsimhits_time.clear();
    fsimhits_rawid.clear();
    fsimhits_energy.clear();
    fsimhits_energy_coded.clear();
    fsimhits_EMfrac.clear();
    fsimhits_bx.clear();
    fsimhits_event.clear();
    fsimhits_ieta.clear();
    fsimhits_iphi.clear();
    fsimhits_eta.clear();
    fsimhits_phi.clear();
    fsimhits_depth.clear();
    fsimhits_subdet.clear();

    double total_energy_of_interest = 0;
    if (!cPrintCSVFile && cPrintAllEvents) cout << "looping on g4SimHits, total " << simhits->size() << ":" << endl;
    for(auto simhit : *simhits)
    {
      // compute
      HcalDetId hcalid( HcalHitRelabeller::relabel( simhit.id(), hcons) );
      int subdet = hcalid.subdet();
      int ietaAbs = hcalid.ietaAbs();
      double EMfrac = simhit.energyEM() / simhit.energy();
      float samplingFactor = 1;
      if(subdet == 1 && ietaAbs-1 < (int)samplingFactors_hb.size()) samplingFactor = samplingFactors_hb.at(ietaAbs-1);
      if(subdet == 2 && ietaAbs-16 < (int)samplingFactors_he.size()) samplingFactor = samplingFactors_he.at(ietaAbs-16);
      std::vector<HcalTrigTowerDetId> ids = tgeo->towerIds(hcalid);  
      int compressed_energy_1 = (outcoder_->compress(ids[0], int(simhit.energy() * samplingFactor), false)).compressedEt();
      double corr_energy = simhit.energy() * samplingFactor;

      // ieta,iphi -> eta,phi
      std::pair<double,double> etaphi = hcons->getEtaPhi(hcalid.subdet(), hcalid.ieta(), hcalid.iphi());
      double eta = etaphi.first;
      double phi = etaphi.second;

      // fill branches
      fsimhits_rawid.push_back(simhit.id());
      fsimhits_time.push_back(simhit.time());
      fsimhits_energy.push_back(simhit.energy() * samplingFactor);
      fsimhits_energy_coded.push_back(compressed_energy_1);
      fsimhits_bx.push_back(simhit.eventId().bunchCrossing());
      fsimhits_event.push_back(simhit.eventId().event());
      fsimhits_ieta.push_back(hcalid.ieta());
      fsimhits_iphi.push_back(hcalid.iphi());
      fsimhits_eta.push_back(eta);
      fsimhits_phi.push_back(phi);
      fsimhits_depth.push_back(hcalid.depth());
      fsimhits_subdet.push_back(hcalid.subdet());
      fsimhits_EMfrac.push_back(EMfrac);

      // fill 2d histos
      if(ids.size() == 1) {
        fsimhit_2d->Fill(ids[0].ieta(), ids[0].iphi(), corr_energy);
      } 
      if(ids.size() == 2) {
        if (hcalid.subdet() == 1) { // HB, split energy
          fsimhit_2d->Fill(ids[0].ieta(), ids[0].iphi(), corr_energy/2.0);
          fsimhit_2d->Fill(ids[1].ieta(), ids[1].iphi(), corr_energy/2.0);
        }
        if (hcalid.subdet() == 2) { // HE, share energy
          fsimhit_2d->Fill(ids[0].ieta(), ids[0].iphi(), corr_energy);
          fsimhit_2d->Fill(ids[1].ieta(), ids[1].iphi(), corr_energy);
        }
        if (hcalid.subdet() == 4) { // HF, undefined
          fsimhit_2d->Fill(ids[0].ieta(), ids[0].iphi(), corr_energy);
        }
      }
      if(ids.size() >= 3) {
        fsimhit_2d->Fill(ids[0].ieta(), ids[0].iphi(), corr_energy);
      }

      // debug
      if (!cPrintCSVFile && (hcalid.ieta() == cTowerIeta && hcalid.iphi() == cTowerIphi)) {
        cout << "found simhit in g4SimHits for " << cTowerIeta << " " << cTowerIphi << endl;
        cout << "energy " << 1000*simhit.energy()*samplingFactor << " MeV, emfrac: " << EMfrac << "\n" << endl;
        total_energy_of_interest += simhit.energy()*samplingFactor;
      }
      if (!cPrintCSVFile && cPrintAllEvents) {
        cout << "" << endl;
        cout << "energy " << simhit.energy() << endl;
        cout << "energyEM " << simhit.energyEM() << endl;
        cout << "id " << simhit.id() << endl;
        cout << "bx " << simhit.eventId().bunchCrossing() << endl;
        cout << "event " << simhit.eventId().event() << endl;
        cout << "time " << simhit.time() << endl;
        cout << "ieta " << hcalid.ieta() << endl;
        cout << "iphi " << hcalid.iphi() << endl;
        cout << "depth " << hcalid.depth() << endl;
        cout << "subdet " << hcalid.subdet() << endl;
      }

      if (cPrintCSVFile) {
        cout << fRunNum << " " << fLumiNum << " " << fEventNum << " | ";
        // name
        if(ids.size() == 1) { cout << "simHit | "; }
        else if(ids.size() == 2) { cout << "simHitDouble | "; }
        else { cout << "simHitOther | "; }
        // energy
        cout << simhit.energy() << " " << samplingFactor << " | ";
        // coded energy
        //int total_samples = 1;
        //int soi = 0;
        //IntegerCaloSamples linear_input(ids[0], total_samples);
        //linear_input.setPresamples(soi);
        //linear_input[soi] = simhit.energy() * samplingFactor;
        //std::vector<int> finegrain_input(total_samples,false); // I think this works, but ulimately need to construct the fg map just like Algo code
        //HcalTriggerPrimitiveDigi compressed_output(ids[0]);
        //outcoder_->compress(linear_input, finegrain_input, compressed_output);
        //int compressed_energy = compressed_output.SOI_compressedEt();
        int compressed_energy = (outcoder_->compress(ids[0], int(simhit.energy() * samplingFactor), false)).compressedEt();
        cout << compressed_energy << " | ";
        // location
        cout << hcalid.ieta() << " " << hcalid.iphi() << " " << hcalid.depth() << " " << hcalid.subdet();
        if(ids.size() == 1) {
          cout << " | " << ids[0].ieta() << " " << ids[0].iphi();
        } else if(ids.size() == 2) {
          cout << " | " << ids[0].ieta() << " " << ids[0].iphi();
          cout << " | " << ids[1].ieta() << " " << ids[1].iphi();
        } else {
          cout << "";
        }
        cout << endl;
      }
    } // end loop on simhits
    if (!cPrintCSVFile && cPrintAllEvents) cout << "done looping on simhits collection" << endl;
    if (!cPrintCSVFile && cTowerIeta != 0 && cTowerIphi != 0) cout << "summed simhit energy in " << cTowerIeta << " " << cTowerIphi << ": " << total_energy_of_interest << " GeV\n" << endl;

    // do compressed simhits 2d histo
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      int bin = fsimhit_2d->FindBin(itr->id().ieta(), itr->id().iphi());
      double linear = fsimhit_2d->GetBinContent(bin);
      double compressed = (outcoder_->compress(itr->id(), int(linear), false)).compressedEt();
      fsimhit_2d_coded->SetBinContent(bin, compressed);
    }

    // dump to CSV summed simhits
    if (cPrintCSVFile) 
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      double simhit = fsimhit_2d->GetBinContent( fsimhit_2d->FindBin(itr->id().ieta(), itr->id().iphi()) );
      double simhit_coded = fsimhit_2d_coded->GetBinContent( fsimhit_2d_coded->FindBin(itr->id().ieta(), itr->id().iphi()) );
      int trig_ieta = itr->id().ieta();
      int trig_iphi = itr->id().iphi();
      cout << fRunNum << " " << fLumiNum << " " << fEventNum << " | ";
      // name
      cout << "simHitSum | ";
      // energy
      cout << simhit << " | ";
      cout << simhit_coded;
      // location
      cout << " | " << trig_ieta << " " << trig_iphi;
      cout << endl;
    }

    // build reverse map
    std::map<std::tuple<int,int,int>, vector<int>> outcoder_revmap;
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      const std::vector<unsigned int> full_lut = outTranscoder_cast->getCompressionLUT(itr->id());
      for(unsigned int i = 0; i < full_lut.size(); i++) {
        int outcoder_from = i;
        int outcoder_to = full_lut[i];
        std::tuple<int,int,int> arg = std::make_tuple(itr->id().ieta(), itr->id().iphi(), outcoder_to);
        if(outcoder_revmap.find( arg ) == outcoder_revmap.end()) {
          outcoder_revmap[arg] = vector<int>();
          outcoder_revmap[arg].push_back(outcoder_from);
        } else {
          outcoder_revmap[arg].push_back(outcoder_from);
        }
      }
      // build reverse map
      /*for(int r = 0; r <= MAX_LINEAR_TP; r++) {
        int outcoder_from = r;
        HcalTriggerPrimitiveSample result = outcoder_->compress(itr->id(), outcoder_from, false);
        int outcoder_to = result.compressedEt();
        std::tuple<int,int,int> arg = std::make_tuple(itr->id().ieta(), itr->id().iphi(), outcoder_to);
        if(outcoder_revmap.find( arg ) == outcoder_revmap.end()) {
          //cout << " new vector" << endl;
          outcoder_revmap[arg] = vector<int>();
          outcoder_revmap[arg].push_back(outcoder_from);
        } else {
          //cout << " add to existting vector" << endl;
          outcoder_revmap[arg].push_back(outcoder_from);
        }
      }*/
    }
    std::map<std::tuple<int,int,int>, double> outcoder_revmap_avg;
    for(auto itr = outcoder_revmap.begin(); itr != outcoder_revmap.end(); itr++) {
      std::tuple<int,int,int> key = itr->first;
      auto vals = itr->second;
      double avg = 0;
      for(auto val = vals.begin(); val != vals.end(); val++) {
        avg += double(*val);
      }
      avg = avg / double(vals.size());
      outcoder_revmap_avg[key] = avg;
    }

    // optionally dump to stdout example of reverse mapping
    bool print_example = false;
    if (print_example) { 
      int compressed = 10;
      int ip = 51;
      cout << "\nexample:" << endl;
      for(int ie = -50; ie < 50; ie++) {
          if(outcoder_revmap.find(std::make_tuple(ie,ip,compressed)) == outcoder_revmap.end())
            continue;
          auto entries = outcoder_revmap.at(std::make_tuple(ie,ip,compressed));
          cout << ie << " : ";
          for(auto loop = entries.begin(); loop != entries.end(); loop++)
          {
            cout << *loop << " ";
          }
          cout << endl;
      }
      cout << "\nafter average:" << endl;
      for(int ie = -50; ie < 50; ie++) {
          if(outcoder_revmap_avg.find(std::make_tuple(ie,ip,compressed)) == outcoder_revmap_avg.end())
            continue;
          cout << ie << " : " << outcoder_revmap_avg.at(std::make_tuple(ie,ip,compressed)) << endl;
      }
      cout << "\nofficial:" << endl;
      for(int ie = -50; ie < 50; ie++) {
          cout << ie << " : " <<  outTranscoder_cast->hcaletValue(ie, ip, 0, compressed) << " " << outTranscoder_cast->hcaletValue(ie, ip, 1, compressed) << endl;
      }
    }

    // big test
    if (false) {
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      int ieta_test = 20;
      int iphi_test = 10;
      int MAX_LINEAR = 1024;
      int MAX_COMPRESSED = 100;
      if(itr->id().ieta() != ieta_test || itr->id().iphi() != iphi_test) continue;
      cout << "working on this tower, ieta " << ieta_test << " iphi " << iphi_test << endl;
   
      cout << "here is the full compression map: " << endl;
      const std::vector<unsigned int> full_lut = outTranscoder_cast->getCompressionLUT(itr->id());
      for(unsigned int i = 0; i < full_lut.size(); i++) {
        cout << "linear " << i << " -> " << full_lut[i] << endl;
      }

      cout << "here, compress each linear integer successively at this tower:" << endl;
      for(int i = 0; i < MAX_LINEAR; i++) {
        cout << "linear " << i << " -> " << (outcoder_->compress(itr->id(), i, false)).compressedEt() << endl;
      }

      cout << "here is what my reverse coder says for this tower:" << endl;
      for(int i = 0; i < MAX_COMPRESSED; i++) {
        cout << "compressed " << i << " -> ";
        if(outcoder_revmap.find(std::make_tuple(ieta_test,iphi_test,i)) == outcoder_revmap.end()) continue;
        auto entries = outcoder_revmap.at(std::make_tuple(ieta_test,iphi_test,i));
        for(auto loop = entries.begin(); loop != entries.end(); loop++) {
          cout << *loop << " ";
        }   
        cout << endl;
      }

      cout << "here is what my reverse coder avg says for this tower:" << endl;
      for(int i = 0; i < MAX_COMPRESSED; i++) {
        if(outcoder_revmap_avg.find(std::make_tuple(ieta_test,iphi_test,i)) == outcoder_revmap_avg.end()) continue;
        cout << "compressed " << i << " -> " << outcoder_revmap_avg.at(std::make_tuple(ieta_test,iphi_test,i)) << endl;
      }

      cout << "here is what the official reverser says for this tower:" << endl;
      for(int i = 0; i < MAX_COMPRESSED; i++) {
        cout << "compressed " << i << " -> " << outTranscoder_cast->hcaletValue(ieta_test, iphi_test, 0, i) << " " << outTranscoder_cast->hcaletValue(ieta_test, iphi_test, 1, i) << endl;
      }
    }
    }

    if (!cPrintCSVFile && cPrintAllEvents) cout << "looping on HcalTrigPrimDigiCollection, total " << emulTPs->size() << ":" << endl;
    // loop on TP collection
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      // compute
      int ieta = itr->id().ieta();
      int iphi = itr->id().iphi();
      int depth = itr->id().depth();
      int subdet = itr->id().subdet();
      int version = itr->id().version();
      uint32_t rawid = itr->id().rawId();
      int soi = itr->presamples();
      int size = itr->size();
      int soi_et = itr->SOI_compressedEt();
      bool soi_finegrain = itr->SOI_fineGrain();
      //cout << "trying to look up ieta, iphi, et " << ieta << " " << iphi << " " << soi_et << endl;
      double soi_et_linear = outcoder_revmap_avg.at(std::make_tuple(ieta, iphi, soi_et));
      //double soi_et_linear = outTranscoder_cast->hcaletValue(ieta, iphi, version, soi_et);

      // ieta,iphi -> eta,phi
      int VERSION = 0; // but maybe need 1?
      std::pair<double,double> etaphi = hcons->getEtaPhi(1, ieta, iphi);
      double eta = etaphi.first;
      double phi = etaphi.second;
      double eta_t_low = 0;
      double eta_t_high = 0;
      tgeo->towerEtaBounds(ieta, VERSION, eta_t_low, eta_t_high);
      double eta_t = (eta_t_high + eta_t_low) / 2.0;
      double phi_seg = 2.0*3.1415926535 / 72;
      double phi_t = (iphi-1) * phi_seg;
      if (ieta >= tgeo->firstHFTower(VERSION)) {
        eta = 0;
        phi = 0;
        eta_t = 0;
        phi_t = 0;
      }

      // fill branches
      fiphi.push_back(iphi);
      fieta.push_back(ieta);
      feta.push_back(eta);
      fphi.push_back(phi);
      feta_t.push_back(eta_t);
      fphi_t.push_back(phi_t);
      fdepth.push_back(depth);
      fsubdet.push_back(subdet);
      fversion.push_back(version);
      frawid.push_back(rawid);
      fsoi.push_back(soi);
      fsize.push_back(size);
      fEt_soi.push_back(soi_et);
      fEt_soi.push_back(soi_et_linear);
      fEt_s0.push_back(itr->sample(0).compressedEt());
      fEt_s1.push_back(itr->sample(1).compressedEt());
      fEt_s2.push_back(itr->sample(2).compressedEt());
      fEt_s3.push_back(itr->sample(3).compressedEt());
      for (int j = 0; j < size; j++) {
        bool j_finegrain = itr->sample(j).fineGrain();
        ffg_all.push_back(j_finegrain);
      }

      // fill histos
      ftp_2d->Fill(ieta, iphi, soi_et);
      ftp_2d_linear->Fill(ieta, iphi, soi_et_linear);

      // debug
      if (!cPrintCSVFile && cPrintAllEvents) {
        cout << "  \nieta, iphi, depth, subdet, rawid: " << ieta << " " << iphi << " " << depth << " " << subdet << " " << rawid << endl;
        cout << "  soi, size: " << soi << " " << size << endl;
        cout << "  soi_Et, soi_finegrain: " << soi_et << " " << soi_finegrain << endl;
        cout << "  all samples:" << endl;
        for (int j = 0; j < size; j++) {
          int j_et = itr->sample(j).compressedEt();
          bool j_finegrain = itr->sample(j).fineGrain();
          uint16_t j_raw = itr->sample(j).raw();
          cout << "    Et, finegrain, raw: " << j_et << " " << j_finegrain << " " << j_raw << endl;
        }
      }
      if (!cPrintCSVFile && ieta == cTowerIeta && iphi == cTowerIphi) {
        cout << "found the output entry of the EDProducer for " << cTowerIeta << " " << cTowerIphi << endl;
        cout << "here are the compressed energies:" << endl;
        for (int j = 0; j < size; j++) {
          int j_et = itr->sample(j).compressedEt();
          cout << j_et << " ";
        }
        cout << endl;  
      }

      if (cPrintCSVFile)
      {
        cout << fRunNum << " " << fLumiNum << " " << fEventNum << " | ";
        cout << "TP | ";
        for (int j = 0; j < size; j++) {
          int j_et = itr->sample(j).compressedEt();
          cout << j_et << " ";
        }
        cout << "| ";
        for (int j = 0; j < size; j++) {
          int j_et = itr->sample(j).compressedEt();
          double linear_e = outcoder_revmap_avg.at(std::make_tuple(ieta, iphi, j_et));
          if(j_et == 0) linear_e = 0;
          cout << linear_e << " ";
        }
        cout << "| " << ieta << " " << iphi;
        cout << endl;
      }
    } // end loop on TPs
    if (!cPrintCSVFile && cPrintAllEvents) cout << "done loop on the TP collection" << endl;
 
    // fill tree
    fTree->Fill();

    // make comparison graphs
    int num_points_count = 0;
    for (HcalTrigPrimDigiCollection::const_iterator itr = emulTPs->begin(); itr != emulTPs->end(); ++itr)
    {
      num_points_count++;
      double tp = ftp_2d->GetBinContent( ftp_2d->FindBin(itr->id().ieta(), itr->id().iphi()) );
      double simhit = fsimhit_2d->GetBinContent( fsimhit_2d->FindBin(itr->id().ieta(), itr->id().iphi()) );

      double tp_linear = ftp_2d_linear->GetBinContent( ftp_2d_linear->FindBin(itr->id().ieta(), itr->id().iphi()) );
      double simhit_coded = fsimhit_2d_coded->GetBinContent( fsimhit_2d_coded->FindBin(itr->id().ieta(), itr->id().iphi()) );

      flinear_compare->SetPoint(num_points_count, simhit, tp_linear);
      fcompressed_compare->SetPoint(num_points_count, simhit_coded, tp);

      flinear_compare_hist->Fill(simhit, tp_linear);
      fcompressed_compare_hist->Fill(simhit_coded, tp);
    }
}

void
L1TPAnalyzer::beginJob()
{
}

void
L1TPAnalyzer::endJob()
{
}

void
L1TPAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TPAnalyzer);
