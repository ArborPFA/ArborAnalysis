/**
 *  @file   PandoraAnalysis/src/CalibrationPfoAnalysis.cc
 * 
 *  @brief  Implementation of the calibration pfo analysis class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCGenericObject.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/CellIDDecoder.h"

#include "CalorimeterHitType.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "CalibrationPfoAnalysis.h"

#include <cmath>
#include <set>

CalibrationPfoAnalysis aCalibrationPfoAnalysis;

//------------------------------------------------------------------------------------------------------------------------------------------

CalibrationPfoAnalysis::CalibrationPfoAnalysis() :
    Processor("CalibrationPfoAnalysis")
{
    _description = "CalibrationPfoAnalysis analyses output of simple clustering for single particle calibration";

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputParticleCollections",
                             "Names of input reconstructed particle collections",
                             m_inputParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                             "MCParticleCollections",
                             "Names of mc particle collections",
                             m_mcParticleCollections,
                             StringVector());

    std::string rootFile("CalibrationPFOAnalysis.root");
    registerProcessorParameter( "RootFile",
                                "Name of the output root file",
                                m_rootFile,
                                rootFile);

    registerProcessorParameter( "Printing",
                                "Set the debug print level",
                                m_printing,
                                (int)0);

    registerProcessorParameter( "UseSemiDigital",
                                "Whether to analyse semi-digital hadronic hits",
                                m_useSemiDigital,
                                (int)1);

    FloatVector semiDigitalThresholds;
    semiDigitalThresholds.push_back(1.f);
    semiDigitalThresholds.push_back(2.f);
    semiDigitalThresholds.push_back(3.f);
    registerProcessorParameter( "SemiDigitalThresholds",
                                "Set the semi-digital threshold to read from hit->getEnergy()",
                                m_semiDigitalThresholds,
                                semiDigitalThresholds);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::init()
{
    m_nRun = 0;
    m_nEvt = 0;

    m_nRunSum = 0;
    m_nEvtSum = 0;

    this->Clear();

    m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

    m_tree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
    m_tree->SetDirectory(m_pTFile);

    m_tree->Branch("run", &m_nRun, "run/I");
    m_tree->Branch("event", &m_nEvt, "event/I");
    m_tree->Branch("ecalHitEnergies", &m_ecalHitEnergies);
    m_tree->Branch("hcalHitEnergies", &m_hcalHitEnergies);
    m_tree->Branch("muonHitEnergies", &m_muonHitEnergies);
    m_tree->Branch("otherHitEnergies", &m_otherHitEnergies);
    m_tree->Branch("nEcalHits", &m_nEcalHits);
    m_tree->Branch("nAhcalHits", &m_nAhcalHits);
    m_tree->Branch("nMuonHits", &m_nMuonHits);
    m_tree->Branch("nOtherHits", &m_nOtherHits);

    if(m_useSemiDigital)
    {
      m_tree->Branch("sdhcalThreshold", &m_sdhcalThreshold);
      m_tree->Branch("hcalOtherHitEnergies", &m_hcalOtherHitEnergies);
      m_tree->Branch("otherEnergySum", &m_otherEnergySum, "otherEnergySum/F");
      m_tree->Branch("hcalOtherEnergySum", &m_hcalOtherEnergySum, "hcalOtherEnergySum/F");
      m_tree->Branch("nSDHcalHits", &m_nSDHcalHits);
      m_tree->Branch("nSDHcalHits1", &m_nSDHcalHits1);
      m_tree->Branch("nSDHcalHits2", &m_nSDHcalHits2);
      m_tree->Branch("nSDHcalHits3", &m_nSDHcalHits3);
    }
    else
    {
      m_tree->Branch("hcalHitEnergies", &m_hcalHitEnergies);
      m_tree->Branch("hcalEnergySum", &m_hcalEnergySum, "hcalEnergySum/F");
    }

    m_tree->Branch("ecalLayerIds", &m_ecalLayerIds);
    m_tree->Branch("hcalLayerIds", &m_hcalLayerIds);
    m_tree->Branch("sdhcalLayerIds", &m_sdhcalLayerIds);
    m_tree->Branch("hcalOtherLayerIds", &m_hcalOtherLayerIds);
    m_tree->Branch("muonLayerIds", &m_muonLayerIds);
    m_tree->Branch("otherLayerIds", &m_otherLayerIds);
    m_tree->Branch("innerEcalLayer", &m_innerEcalLayer);
    m_tree->Branch("outerEcalLayer", &m_outerEcalLayer);
    m_tree->Branch("innerHcalLayer", &m_innerHcalLayer);
    m_tree->Branch("outerHcalLayer", &m_outerHcalLayer);
    m_tree->Branch("innerMuonLayer", &m_innerMuonLayer);
    m_tree->Branch("outerMuonLayer", &m_outerMuonLayer);
    m_tree->Branch("innerOtherLayer", &m_innerOtherLayer);
    m_tree->Branch("outerOtherLayer", &m_outerOtherLayer);
    m_tree->Branch("muonEnergySum", &m_muonEnergySum, "muonEnergySum/F");
    m_tree->Branch("ecalEnergySum", &m_ecalEnergySum, "ecalEnergySum/F");
    m_tree->Branch("pfoTotalEnergy", &m_pfoTotalEnergy, "mpfoTotalEnergy/F");
    m_tree->Branch("mcNParticles", &m_mcNParticles, "mcNParticles/I");
    m_tree->Branch("mcNPhotons", &m_mcNPhotons, "mcNPhotons/I");
    m_tree->Branch("mcNTracks", &m_mcNTracks, "mcNTracks/I");
    m_tree->Branch("mcNNeutralHadrons", &m_mcNNeutralHadrons, "mcNNeutralHadrons/I");
    m_tree->Branch("mcTotalEnergy", &m_mcTotalEnergy, "mcTotalEnergy/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::processRunHeader(lcio::LCRunHeader */*pLCRunHeader*/)
{
    m_nRun = 0;
    m_nEvt = 0;

    ++m_nRunSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::processEvent(EVENT::LCEvent *pLCEvent)
{
    m_nRun = pLCEvent->getRunNumber();
    m_nEvt = pLCEvent->getEventNumber();
    ++m_nEvtSum;

    if ((m_nEvtSum % 100) == 0)
        std::cout << " processed events: " << m_nEvtSum << std::endl;

    this->Clear();
    this->ExtractCollections(pLCEvent);
    this->PerformPfoAnalysis();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::check(EVENT::LCEvent */*pLCEvent*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::end()
{
    if (m_printing > -1)
    {
        std::cout << "PfoAnalysis::end() " << this->name() << " processed " << m_nEvtSum << " events in " << m_nRunSum << " runs " << std::endl
                  << "Rootfile: " << m_rootFile.c_str() << std::endl;
    }

    m_tree->Write();

    m_pTFile->Write();
    m_pTFile->Close();
    delete m_pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::Clear()
{
    m_pfoVector.clear();
    m_pfoTargetVector.clear();

    m_ecalHitEnergies.clear();
    m_hcalHitEnergies.clear();
    m_muonHitEnergies.clear();
    m_hcalOtherHitEnergies.clear();
    m_otherHitEnergies.clear();
    m_sdhcalThreshold.clear();
    m_ecalLayerIds.clear();
    m_hcalLayerIds.clear();
    m_sdhcalLayerIds.clear();
    m_hcalOtherLayerIds.clear();
    m_muonLayerIds.clear();
    m_otherLayerIds.clear();

    m_nEcalHits = 0;
    m_nAhcalHits = 0;
    m_nSDHcalHits = 0;
    m_nSDHcalHits1 = 0;
    m_nSDHcalHits2 = 0;
    m_nSDHcalHits3 = 0;
    m_nMuonHits = 0;
    m_nOtherHits = 0;

    m_innerEcalLayer = 100;
    m_outerEcalLayer = -1;
    m_innerHcalLayer = 100;
    m_outerHcalLayer = -1;
    m_innerMuonLayer = 100;
    m_outerMuonLayer = -1;
    m_innerOtherLayer = 100;
    m_outerOtherLayer = -1;

    m_ecalEnergySum = 0.f;
    m_hcalEnergySum = 0.f;
    m_hcalOtherEnergySum = 0.f;
    m_muonEnergySum = 0.f;
    m_otherEnergySum = 0.f;
    m_pfoTotalEnergy = 0.f;

    m_mcNParticles = 0;
    m_mcNPhotons = 0;
    m_mcNTracks = 0;
    m_mcNNeutralHadrons = 0;

    m_mcTotalEnergy = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
  // Extract mc particle collection
  std::set<MCParticle*> pfoTargetList;

  for (StringVector::const_iterator iter = m_mcParticleCollections.begin(), iterEnd = m_mcParticleCollections.end();
      iter != iterEnd; ++iter)
  {
      try
      {
          const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

          for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
          {
              MCParticle *pMCParticle = dynamic_cast<MCParticle*>(pLCCollection->getElementAt(i));

              if (NULL == pMCParticle)
                  throw EVENT::Exception("Collection type mismatch");

              const float innerRadius(std::sqrt(pMCParticle->getVertex()[0] * pMCParticle->getVertex()[0] + pMCParticle->getVertex()[1] * pMCParticle->getVertex()[1] + pMCParticle->getVertex()[2] * pMCParticle->getVertex()[2]));
              const float outerRadius(std::sqrt(pMCParticle->getEndpoint()[0] * pMCParticle->getEndpoint()[0] + pMCParticle->getEndpoint()[1] * pMCParticle->getEndpoint()[1] + pMCParticle->getEndpoint()[2] * pMCParticle->getEndpoint()[2]));
              const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));

              if ((pfoTargetList.find(pMCParticle) == pfoTargetList.end()) && (outerRadius > 500.f) && (innerRadius <= 500.f) &&
                  (momentum > 0.01f) && !((pMCParticle->getPDG() == 2212 || pMCParticle->getPDG() == 2112) && (pMCParticle->getEnergy() < 1.2f)))
              {
                  pfoTargetList.insert(pMCParticle);
                  m_pfoTargetVector.push_back(pMCParticle);
              }
          }
      }
      catch (...)
      {
          streamlog_out(WARNING) << "Could not extract mc particle information" << std::endl;
      }
  }

    // Extract reconstructed particle collection
    for (StringVector::const_iterator iter = m_inputParticleCollections.begin(), iterEnd = m_inputParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));

                if (NULL == pReconstructedParticle)
                    throw EVENT::Exception("Collection type mismatch");

                m_pfoVector.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input particle collection: " << *iter << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationPfoAnalysis::PerformPfoAnalysis()
{
    CellIDDecoder<CalorimeterHit> cellIDDecoder("M:3,S-1:3,I:9,J:9,K-1:6");

    // Extract quantities relating to reconstructed pfos
    for (ParticleVector::const_iterator iter = m_pfoVector.begin(), iterEnd = m_pfoVector.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pPfo = *iter;

        m_pfoTotalEnergy += pPfo->getEnergy();

        const ClusterVec &clusterVec = pPfo->getClusters();

        for (ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
        {
            const CalorimeterHitVec &calorimeterHitVec = (*iter)->getCalorimeterHits();

            for (CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
            {
                EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                int layer = cellIDDecoder(pCalorimeterHit)["K-1"];
                const float hitEnergy(pCalorimeterHit->getEnergy());
                CHT cht(pCalorimeterHit->getType());

                if (cht.is(CHT::ecal))
                {
                  if(m_innerEcalLayer > layer)
                    m_innerEcalLayer = layer;

                  if(m_outerEcalLayer < layer)
                    m_outerEcalLayer = layer;

                  m_nEcalHits++;
                  m_ecalHitEnergies.push_back(hitEnergy);
                  m_ecalLayerIds.push_back(layer);
                  m_ecalEnergySum += hitEnergy;
                }
                else if (cht.is(CHT::hcal))
                {
                  if(m_innerHcalLayer > layer)
                    m_innerHcalLayer = layer;

                  if(m_outerHcalLayer < layer)
                    m_outerHcalLayer = layer;

                  if(m_useSemiDigital)
                  {
                    if(m_semiDigitalThresholds.at(0) == hitEnergy)
                    {
                      m_nSDHcalHits++;
                      m_nSDHcalHits1++;
                      m_sdhcalThreshold.push_back(1);
                      m_sdhcalLayerIds.push_back(layer);
                    }
                    else if(m_semiDigitalThresholds.at(1) == hitEnergy)
                    {
                      m_nSDHcalHits++;
                      m_nSDHcalHits2++;
                      m_sdhcalThreshold.push_back(2);
                      m_sdhcalLayerIds.push_back(layer);
                    }
                    else if(m_semiDigitalThresholds.at(2) == hitEnergy)
                    {
                      m_nSDHcalHits++;
                      m_nSDHcalHits3++;
                      m_sdhcalThreshold.push_back(3);
                      m_sdhcalLayerIds.push_back(layer);
                    }
                    else
                    {
                      m_nAhcalHits++;
                      m_hcalOtherHitEnergies.push_back(hitEnergy);
                      m_hcalOtherEnergySum += hitEnergy;
                      m_hcalOtherLayerIds.push_back(layer);
                    }
                  }
                  else
                  {
                    m_hcalHitEnergies.push_back(hitEnergy);
                    m_hcalLayerIds.push_back(layer);
                    m_nAhcalHits++;
                  }
                }
                else if (cht.is(CHT::muon))
                {
                  if(m_innerMuonLayer > layer)
                    m_innerMuonLayer = layer;

                  if(m_outerMuonLayer < layer)
                    m_outerMuonLayer = layer;

                  m_muonHitEnergies.push_back(hitEnergy);
                  m_muonEnergySum += hitEnergy;
                  m_muonLayerIds.push_back(layer);
                  m_nMuonHits++;
                }
                else
                {
                  if(m_innerOtherLayer > layer)
                    m_innerOtherLayer = layer;

                  if(m_outerOtherLayer < layer)
                    m_outerOtherLayer = layer;

                  m_otherHitEnergies.push_back(hitEnergy);
                  m_otherEnergySum += hitEnergy;
                  m_otherLayerIds.push_back(layer);
                  m_nOtherHits++;
                }
            }
        }
    }

    // Perform checks on pfo targets : pre-interaction
    for (MCParticleVector::const_iterator iter = m_pfoTargetVector.begin(), iterEnd = m_pfoTargetVector.end(); iter != iterEnd; ++iter)
    {
        MCParticle *pMCParticle = *iter;

        ++m_mcNParticles;
        m_mcTotalEnergy += pMCParticle->getEnergy();

        if (22 == pMCParticle->getPDG())
        {
            ++m_mcNPhotons;
        }
        else if ((11 == std::abs(pMCParticle->getPDG())) || (13 == std::abs(pMCParticle->getPDG())) || (211 == std::abs(pMCParticle->getPDG()))) // TODO, more options here?
        {
            ++m_mcNTracks;
        }
        else
        {
            ++m_mcNNeutralHadrons;
        }
    }

    m_tree->Fill();
}
