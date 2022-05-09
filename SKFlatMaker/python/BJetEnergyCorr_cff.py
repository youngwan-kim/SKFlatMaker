from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
jetCorrFactors = patJetCorrFactors.clone(
    src='slimmedJets',
    levels = cms.vstring(
        'L1FastJet',
        'L2Relative',
        'L3Absolute',
        'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *
updatedJets = updatedPatJets.clone(
    addBTagInfo = False,
    jetSource = 'slimmedJets',
    jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactors") ),
)

bJetVars = cms.EDProducer(
    "JetRegressionVarProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("updatedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
)

updatedJetsWithUserData = cms.EDProducer(
    "PATJetUserDataEmbedder",
    src = cms.InputTag("updatedJets"),
    userFloats = cms.PSet(
        leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
        leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
        leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
        leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
        leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
        leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
        leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
        leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
        leptonPt = cms.InputTag("bJetVars:leptonPt"),
        vtxPt = cms.InputTag("bJetVars:vtxPt"),
        vtxMass = cms.InputTag("bJetVars:vtxMass"),
        vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
        vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
        ptD = cms.InputTag("bJetVars:ptD"),
        genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
    ),
    userInts = cms.PSet(
        vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
        leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
    ),
)

jetUpdaterSeq = cms.Sequence(
    jetCorrFactors
    + updatedJets
    + bJetVars
    + updatedJetsWithUserData
)

bJetNN = cms.EDProducer(
    "BJetEnergyRegressionMVA",
    backend = cms.string("TF"),
    src = cms.InputTag("updatedJetsWithUserData"),
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    rhosrc = cms.InputTag("fixedGridRhoFastjetAll"),

    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2018.pb"),
    name = cms.string("JetRegNN"),
    isClassifier = cms.bool(False),
    variablesOrder = cms.vstring(["Jet_pt","Jet_eta","rho","Jet_mt","Jet_leadTrackPt","Jet_leptonPtRel","Jet_leptonDeltaR","Jet_neHEF",
                                  "Jet_neEmEF","Jet_vtxPt","Jet_vtxMass","Jet_vtx3dL","Jet_vtxNtrk","Jet_vtx3deL",
                                  "Jet_numDaughters_pt03","Jet_energyRing_dR0_em_Jet_rawEnergy","Jet_energyRing_dR1_em_Jet_rawEnergy",
                                  "Jet_energyRing_dR2_em_Jet_rawEnergy","Jet_energyRing_dR3_em_Jet_rawEnergy","Jet_energyRing_dR4_em_Jet_rawEnergy",
                                  "Jet_energyRing_dR0_neut_Jet_rawEnergy","Jet_energyRing_dR1_neut_Jet_rawEnergy","Jet_energyRing_dR2_neut_Jet_rawEnergy",
                                  "Jet_energyRing_dR3_neut_Jet_rawEnergy","Jet_energyRing_dR4_neut_Jet_rawEnergy","Jet_energyRing_dR0_ch_Jet_rawEnergy",
                                  "Jet_energyRing_dR1_ch_Jet_rawEnergy","Jet_energyRing_dR2_ch_Jet_rawEnergy","Jet_energyRing_dR3_ch_Jet_rawEnergy",
                                  "Jet_energyRing_dR4_ch_Jet_rawEnergy","Jet_energyRing_dR0_mu_Jet_rawEnergy","Jet_energyRing_dR1_mu_Jet_rawEnergy",
                                  "Jet_energyRing_dR2_mu_Jet_rawEnergy","Jet_energyRing_dR3_mu_Jet_rawEnergy","Jet_energyRing_dR4_mu_Jet_rawEnergy",
                                  "Jet_chHEF","Jet_chEmEF","Jet_leptonPtRelInv","isEle","isMu","isOther","Jet_mass","Jet_ptd"]),
    variables = cms.PSet(
        Jet_pt = cms.string("pt*jecFactor('Uncorrected')"),
        Jet_mt = cms.string("mt*jecFactor('Uncorrected')"),
        Jet_eta = cms.string("eta"),
        Jet_mass = cms.string("mass*jecFactor('Uncorrected')"),
        Jet_ptd = cms.string("userFloat('ptD')"),
        Jet_leadTrackPt = cms.string("userFloat('leadTrackPt')"),
        Jet_vtxNtrk = cms.string("userInt('vtxNtrk')"),
        Jet_vtxMass = cms.string("userFloat('vtxMass')"),
        Jet_vtx3dL = cms.string("userFloat('vtx3dL')"),
        Jet_vtx3deL = cms.string("userFloat('vtx3deL')"),
        Jet_vtxPt = cms.string("userFloat('vtxPt')"),
        Jet_leptonPtRel = cms.string("userFloat('leptonPtRelv0')"),
        Jet_leptonPtRelInv = cms.string("userFloat('leptonPtRelInvv0')*jecFactor('Uncorrected')"),
        Jet_leptonDeltaR = cms.string("userFloat('leptonDeltaR')"),
        Jet_neHEF = cms.string("neutralHadronEnergyFraction()"),
        Jet_neEmEF = cms.string("neutralEmEnergyFraction()"),
        Jet_chHEF = cms.string("chargedHadronEnergyFraction()"),
        Jet_chEmEF = cms.string("chargedEmEnergyFraction()"),
        isMu = cms.string("?abs(userInt('leptonPdgId'))==13?1:0"),
        isEle = cms.string("?abs(userInt('leptonPdgId'))==11?1:0"),
        isOther = cms.string("?userInt('leptonPdgId')==0?1:0"),
    ),
    inputTensorName = cms.string("ffwd_inp"),
    outputTensorName = cms.string("ffwd_out/BiasAdd"),
    outputNames = cms.vstring(["bJetNNCorr","bJetNNRes"]),
    outputFormulas = cms.vstring(["at(0)*0.27912887930870056+1.0545977354049683","0.5*(at(2)-at(1))*0.27912887930870056"]),
    nThreads = cms.uint32(1),
    singleThreadPool = cms.string("no_threads"),
)#bJetNN

cJetNN= cms.EDProducer(
    "BJetEnergyRegressionMVA",
    backend = cms.string("TF"),
    src = cms.InputTag("updatedJetsWithUserData"),
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    rhosrc = cms.InputTag("fixedGridRhoFastjetAll"),
    
    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2018.pb"),
name = cms.string("JetRegNN"),
    isClassifier = cms.bool(False),
    variablesOrder = cms.vstring(["Jet_pt","Jet_eta","rho","Jet_mt","Jet_leadTrackPt","Jet_leptonPtRel","Jet_leptonDeltaR",
                                  "Jet_neHEF","Jet_neEmEF","Jet_vtxPt","Jet_vtxMass","Jet_vtx3dL","Jet_vtxNtrk","Jet_vtx3deL",
                                  "Jet_numDaughters_pt03","Jet_chEmEF","Jet_chHEF", "Jet_ptd","Jet_mass",
                                  "Jet_energyRing_dR0_em_Jet_rawEnergy","Jet_energyRing_dR1_em_Jet_rawEnergy",
                                  "Jet_energyRing_dR2_em_Jet_rawEnergy","Jet_energyRing_dR3_em_Jet_rawEnergy","Jet_energyRing_dR4_em_Jet_rawEnergy",
                                  "Jet_energyRing_dR0_neut_Jet_rawEnergy","Jet_energyRing_dR1_neut_Jet_rawEnergy","Jet_energyRing_dR2_neut_Jet_rawEnergy",
                                  "Jet_energyRing_dR3_neut_Jet_rawEnergy","Jet_energyRing_dR4_neut_Jet_rawEnergy","Jet_energyRing_dR0_ch_Jet_rawEnergy",
                                  "Jet_energyRing_dR1_ch_Jet_rawEnergy","Jet_energyRing_dR2_ch_Jet_rawEnergy","Jet_energyRing_dR3_ch_Jet_rawEnergy",
                                  "Jet_energyRing_dR4_ch_Jet_rawEnergy","Jet_energyRing_dR0_mu_Jet_rawEnergy","Jet_energyRing_dR1_mu_Jet_rawEnergy",
                                  "Jet_energyRing_dR2_mu_Jet_rawEnergy","Jet_energyRing_dR3_mu_Jet_rawEnergy","Jet_energyRing_dR4_mu_Jet_rawEnergy"]),
    variables = cms.PSet(
        Jet_pt = cms.string("pt*jecFactor('Uncorrected')"),
        Jet_mt = cms.string("mt*jecFactor('Uncorrected')"),
        Jet_eta = cms.string("eta"),
        Jet_mass = cms.string("mass*jecFactor('Uncorrected')"),
        Jet_ptd = cms.string("userFloat('ptD')"),
        Jet_leadTrackPt = cms.string("userFloat('leadTrackPt')"),
        Jet_vtxNtrk = cms.string("userInt('vtxNtrk')"),
        Jet_vtxMass = cms.string("userFloat('vtxMass')"),
        Jet_vtx3dL = cms.string("userFloat('vtx3dL')"),
        Jet_vtx3deL = cms.string("userFloat('vtx3deL')"),
        Jet_vtxPt = cms.string("userFloat('vtxPt')"),
        Jet_leptonPtRel = cms.string("userFloat('leptonPtRelv0')"),
        Jet_leptonPtRelInv = cms.string("userFloat('leptonPtRelInvv0')*jecFactor('Uncorrected')"),
        Jet_leptonDeltaR = cms.string("userFloat('leptonDeltaR')"),
        Jet_neHEF = cms.string("neutralHadronEnergyFraction()"),
        Jet_neEmEF = cms.string("neutralEmEnergyFraction()"),
        Jet_chHEF = cms.string("chargedHadronEnergyFraction()"),
        Jet_chEmEF = cms.string("chargedEmEnergyFraction()"),
        isMu = cms.string("?abs(userInt('leptonPdgId'))==13?1:0"),
        isEle = cms.string("?abs(userInt('leptonPdgId'))==11?1:0"),
        isOther = cms.string("?userInt('leptonPdgId')==0?1:0"),
    ),
    inputTensorName = cms.string("ffwd_inp"),
    outputTensorName = cms.string("ffwd_out/BiasAdd"),
    outputNames = cms.vstring(["cJetNNCorr","cJetNNRes"]),
    outputFormulas = cms.vstring(["at(0)*0.24325256049633026+0.993854820728302","0.5*(at(2)-at(1))*0.24325256049633026"]),
    nThreads = cms.uint32(1),
    singleThreadPool = cms.string("no_threads"),
)#cjetNN

correctionSeq = cms.Sequence(
    bJetNN
    + cJetNN
)

bJetCorrSeq = cms.Sequence(
    jetUpdaterSeq
    + correctionSeq
)
