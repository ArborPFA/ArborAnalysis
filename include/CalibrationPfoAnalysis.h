/**
 *  @file   ArborAnalysis/include/CalibrationPfoAnalysis.h
 * 
 *  @brief  Header file for the calibration pfo analysis class.
 * 
 *  $Log: $
 */

#ifndef CALIBRATION_PFO_ANALYSIS_H
#define CALIBRATION_PFO_ANALYSIS_H 1

#include "EVENT/ReconstructedParticle.h"

#include "marlin/Processor.h"

class TFile;
class TH1F;
class TTree;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CalibrationPfoAnalysis class
 */
class CalibrationPfoAnalysis : public marlin::Processor
{
public:
    /**
     *  @brief  Default constructor
     */
    CalibrationPfoAnalysis();

    /**
     *  @brief  Create new processor
     */
    virtual Processor *newProcessor();

    /**
     *  @brief  Initialize, called at startup
     */
    virtual void init();

    /**
     *  @brief  Process run header
     *
     *  @param  pLCRunHeader the lc run header
     */
    virtual void processRunHeader(lcio::LCRunHeader *pLCRunHeader);

    /**
     *  @brief  Process event, main entry point
     *
     *  @param  pLCEvent the lc event
     */
    virtual void processEvent(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  End, called at shutdown
     */
    virtual void end();

private:
    /**
     *  @brief  Clear current event details
     */
    void Clear();

    /**
     *  @brief  Extract lcio collections
     * 
     *  @param  pLCEvent the lc event
     */
    void ExtractCollections(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Perform pfo analysis
     */
    void PerformPfoAnalysis();

    typedef std::vector<ReconstructedParticle *> ParticleVector;
    typedef std::vector<MCParticle*> MCParticleVector;
    typedef std::vector<std::string> StringVector;
    typedef std::vector<float> FloatVector;
    typedef std::vector<int> IntVector;

    int                 m_nRun;                                 ///< 
    int                 m_nEvt;                                 ///< 

    int                 m_nRunSum;                              ///<
    int                 m_nEvtSum;                              ///<

    StringVector        m_mcParticleCollections;                ///<
    StringVector        m_inputParticleCollections;             ///< 

    int                 m_printing;                             ///<
    int                 m_useSemiDigital;                       ///<
    std::string         m_rootFile;                             ///< 

    FloatVector         m_semiDigitalThresholds;                ///<

    ParticleVector      m_pfoVector;                            ///< 
    MCParticleVector    m_pfoTargetVector;                      ///<


    FloatVector         m_ecalHitEnergies;
    FloatVector         m_hcalHitEnergies;  // only for ahcal
    FloatVector         m_muonHitEnergies;
    FloatVector         m_hcalOtherHitEnergies; // only for sdhcal
    FloatVector         m_otherHitEnergies;
    IntVector           m_sdhcalThreshold;  // only for sdhcal
    IntVector           m_ecalLayerIds;
    IntVector           m_hcalLayerIds;
    IntVector           m_sdhcalLayerIds;
    IntVector           m_hcalOtherLayerIds;
    IntVector           m_muonLayerIds;
    IntVector           m_otherLayerIds;

    int                 m_nEcalHits;
    int                 m_nAhcalHits;
    int                 m_nSDHcalHits;
    int                 m_nSDHcalHits1;
    int                 m_nSDHcalHits2;
    int                 m_nSDHcalHits3;
    int                 m_nMuonHits;
    int                 m_nOtherHits;

    int                 m_innerEcalLayer;
    int                 m_outerEcalLayer;
    int                 m_innerHcalLayer;
    int                 m_outerHcalLayer;
    int                 m_innerMuonLayer;
    int                 m_outerMuonLayer;
    int                 m_innerOtherLayer;
    int                 m_outerOtherLayer;

    float               m_ecalEnergySum;                        ///<
    float               m_hcalEnergySum;                        ///<
    float               m_hcalOtherEnergySum;                   ///<
    float               m_muonEnergySum;                        ///<
    float               m_otherEnergySum;                       ///<
    float               m_pfoTotalEnergy;                       ///<

    int                 m_mcNParticles;                         ///<
    int                 m_mcNPhotons;                           ///<
    int                 m_mcNTracks;                            ///<
    int                 m_mcNNeutralHadrons;                    ///<

    float               m_mcTotalEnergy;                        ///<

    TFile              *m_pTFile;                               ///< 
    TTree              *m_tree;                                 ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *CalibrationPfoAnalysis::newProcessor()
{
    return new CalibrationPfoAnalysis;
}

#endif // #ifndef CALIBRATION_PFO_ANALYSIS_H
