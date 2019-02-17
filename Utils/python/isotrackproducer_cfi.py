import FWCore.ParameterSet.Config as cms

isotrackproducer = cms.EDProducer("IsoTrackProducer",
                                    selectedTracks = cms.InputTag("generalTracks"),
                                    #selectedTracks = cms.InputTag("isolatedTracks"),
                                    selectedElectrons = cms.InputTag("gedGsfElectrons"),
                                    selectedMuons = cms.InputTag("muons"),
                                    selectedPFJets = cms.InputTag("ak4PFJetsCHS"),
                                    selectedPFCand = cms.InputTag("particleFlow"),
                                    minTrackPt = cms.double(15),
                                    maxTrackEta = cms.double(2.4),
                                    coneRelIsoDR = cms.double(0.3),
                                    conePtSumMaxPtFraction = cms.double(0.05),
                                    minTrackJetDR = cms.double(0.5),
                                    minTrackLeptonDR = cms.double(0.15),
                                    RequireNumberOfValidPixelHits = cms.int32(1),
                                    RequireNumberOfValidTrackerHits = cms.int32(7),
                                    maxDxy = cms.double(0.02),
                                    maxDz = cms.double(0.5),
                                    selectedEcalRecHitsEB = cms.InputTag("reducedEcalRecHitsEB"),
                                    selectedEcalRecHitsEE = cms.InputTag("reducedEcalRecHitsEE"),
                                    selectedEcalRecHitsES = cms.InputTag("reducedEcalRecHitsES"),
                                    selectedHcalRecHits = cms.InputTag("reducedHcalRecHits"),
                                    selectedCaloJets = cms.InputTag("ak4CaloJets"),
                                    useCaloJetsInsteadOfHits = cms.bool(False),
                                    minMissingOuterHits = cms.int32(3),
                                    caloEnergyDepositionMaxDR = cms.double(0.5),
                                    caloEnergyDepositionMaxE = cms.double(10),
                                    deadNoisyDR = cms.double(0.05),
                                    #dEdxEstimator = cms.string('dedxHarmonic2'), #original pre-pat
    #from PAT:                                
    dEdxDataStrip = cms.InputTag("dedxHarmonic2"),
    dEdxDataPixel = cms.InputTag("dedxPixelHarmonic2"),
    addPrescaledDeDxTracks = cms.bool(False),
    dEdxHitInfo = cms.InputTag("dedxHitInfo"),
	usePrecomputedDeDxStrip = cms.bool(True),        # if these are set to True, will get estimated DeDx from DeDxData branches    
    usePrecomputedDeDxPixel = cms.bool(True),        # if set to False, will manually compute using dEdxHitInfo	
                                        
                                    doDeDx = cms.bool(True),
                                    PrimaryVertex = cms.InputTag('offlinePrimaryVertices'),
                                    maxChargedPFCandSumDR = cms.double(0.01),
                                    maxNeutralPFCandSumDR = cms.double(0.05),
                                    )


''' from /cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_11/src/PhysicsTools/PatAlgos/python/slimming/isolatedTracks_cfi.py
    tkAssocParamBlock,
    packedPFCandidates = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    generalTracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    caloJets = cms.InputTag("ak4CaloJets"),
    dEdxDataStrip = cms.InputTag("dedxHarmonic2"),
    dEdxDataPixel = cms.InputTag("dedxPixelHarmonic2"),
    dEdxHitInfo = cms.InputTag("dedxHitInfo"),
    dEdxHitInfoPrescale = cms.InputTag("dedxHitInfo","prescale"), 
    addPrescaledDeDxTracks = cms.bool(False),
    usePrecomputedDeDxStrip = cms.bool(True),        # if these are set to True, will get estimated DeDx from DeDxData branches
    usePrecomputedDeDxPixel = cms.bool(True),        # if set to False, will manually compute using dEdxHitInfo
    pT_cut = cms.double(5.0),         # save tracks above this pt
    pT_cut_noIso = cms.double(20.0),  # for tracks with at least this pT, don't apply any iso cut
    pfIsolation_DR = cms.double(0.3),
    pfIsolation_DZ = cms.double(0.1),
    miniIsoParams = cms.vdouble(0.05, 0.2, 10.0), # (minDR, maxDR, kT)
                                                  # dR for miniiso is max(minDR, min(maxDR, kT/pT))
    absIso_cut = cms.double(5.0),
    relIso_cut = cms.double(0.2),
    miniRelIso_cut = cms.double(0.2),

    caloJet_DR = cms.double(0.3),

    pflepoverlap_DR = cms.double(0.001),
    pflepoverlap_pTmin = cms.double(5.0),

    pcRefNearest_DR = cms.double(0.3),
    pcRefNearest_pTmin = cms.double(5.0),

    pfneutralsum_DR = cms.double(0.05),

    saveDeDxHitInfo = cms.bool(True),
    saveDeDxHitInfoCut = cms.string("(%s) || (%s)" % (_susySoftDisappearingTrackCut,_exoHighPtTrackCut)), 
    '''