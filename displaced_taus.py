import ROOT
import random
import numpy as np
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from scipy.constants import c as speed_of_light

# from CMGTools.HNL.samples.signals_2018 import all_signals_t
# 
# allfiles = []
# for isample in all_signals_t:
#     if isample.mass < 4: continue
#     allfiles += isample.files

def isAncestor(a, p):
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

def finalDaughters(gen, daughters=None):
    if daughters is None:
        daughters = []
    for i in range(gen.numberOfDaughters()):
        daughter = gen.daughter(i)
        if daughter.numberOfDaughters() == 0:
            daughters.append(daughter)
        else:
            finalDaughters(daughter, daughters)
    return daughters

branches = [
    'run'      ,  
    'lumi'     , 
    'event'    ,

    'gen_tau_pt'   ,
    'gen_tau_eta'  ,
    'gen_tau_phi'  ,
    'gen_tau_pdgid',
    'gen_tau_q'    ,
    'gen_tau_dm'   ,

    'gen_tau_vis_pt'   ,
    'gen_tau_vis_eta'  ,
    'gen_tau_vis_phi'  ,
    'gen_tau_vis_mass' ,

    'reco_tau_pt'   ,
    'reco_tau_eta'  ,
    'reco_tau_phi'  ,
    'reco_tau_mass' ,
    'reco_tau_pdgid',
    'reco_tau_q'    ,
    'reco_tau_dxy'  ,
    'reco_tau_dz'   ,
    'reco_tau_dm'   ,
#     'reco_tau_dtj'  , # DeepTau doesn't exist in these samples!
#     'reco_tau_dte'  , # DeepTau doesn't exist in these samples!
#     'reco_tau_dtm'  , # DeepTau doesn't exist in these samples!
    'reco_tau_idj'  ,
    'reco_tau_ide'  ,
    'reco_tau_idm'  ,

    'reco_tau_bis_pt'   ,
    'reco_tau_bis_eta'  ,
    'reco_tau_bis_phi'  ,
    'reco_tau_bis_mass' ,
    'reco_tau_bis_pdgid',
    'reco_tau_bis_q'    ,
    'reco_tau_bis_dxy'  ,
    'reco_tau_bis_dz'   ,
    'reco_tau_bis_dm'   ,
#     'reco_tau_bis_dtj'  , # DeepTau doesn't exist in these samples!
#     'reco_tau_bis_dte'  , # DeepTau doesn't exist in these samples!
#     'reco_tau_bis_dtm'  , # DeepTau doesn't exist in these samples!
    'reco_tau_bis_idj'  ,
    'reco_tau_bis_ide'  ,
    'reco_tau_bis_idm'  ,

    'Lxy' , # 2D transverse displacement
    'Lxyz', # 3D displacement
]

handles = OrderedDict()
handles['genp_pruned'] = ('prunedGenParticles' , Handle('std::vector<reco::GenParticle>'))
handles['lhe']         = ('externalLHEProducer', Handle('LHEEventProduct'))
handles['taus']        = ('slimmedTaus'        , Handle('std::vector<pat::Tau>'))
handles['taus_bis']    = ('selectedPatTaus'    , Handle('std::vector<pat::Tau>'))

# output file and tree gymnastics
outfile = ROOT.TFile.Open('displaced_taus_all_in_once.root', 'recreate')
ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill  = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

# random.shuffle(allfiles) #
events = Events(['HNL_miniAOD_rerunTauRECO.root'])

for i, event in enumerate(events):

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    if i%1000==0:
        percentage = float(i)/events.size()*100.
        print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

    # all gen particles
    event.genp = [ip for ip in event.genp_pruned]

    # get the heavy neutrino
    the_hns = [ip for ip in event.genp_pruned if abs(ip.pdgId()) in [9900012, 9990012] and ip.isLastCopy()] # 9900012 is Majorana, 9990012 is Dirac. Dirac comes in two species, particle and anti-particle!
    event.the_hn = the_hns[0] # one per event

    # prompt lepton
    prompt_leptons = []
    prompt_leptons += [ip for ip in event.genp_pruned if abs(ip.pdgId()) in [11,13] and ip.isPromptFinalState() and not isAncestor(event.the_hn, ip)]
    prompt_leptons += [ip for ip in event.genp_pruned if abs(ip.pdgId()) in [15]    and ip.isPromptDecayed()    and not isAncestor(event.the_hn, ip)]
    event.the_pl = prompt_leptons[0]
         
    the_taus = [ip for ip in event.genp_pruned if abs(ip.pdgId())==15 and ip.isLastCopy()]
    
    #print '\tthe taus'
    for itau in the_taus:
        if itau == event.the_pl: continue
        if itau.pt() < 15 or abs(itau.eta())>2.4: continue 
        # print itau.pdgId(), itau.pt(), itau.eta(), itau.phi(),
        daughters_pdgids = [abs(itau.daughter(idau).pdgId()) for idau in xrange(itau.numberOfDaughters())]
        if 11 in daughters_pdgids or 13 in daughters_pdgids:
            continue
        
        # visible gen tau p4
        itau_vis = ROOT.Math.LorentzVector('<ROOT::Math::PtEtaPhiM4D<double> >')()
        for idau in xrange(itau.numberOfDaughters()):
            if abs(itau.daughter(idau).pdgId()) not in [12, 14, 16]:
                itau_vis += itau.daughter(idau).p4()
    
        # gen decay mode
        itau.gen_dm = tauDecayModes.genDecayModeInt([d for d in finalDaughters(itau) if abs(d.pdgId()) not in [12, 14, 16]])

        # reco tau
        reco_tau = None
        match = None
        dr2 = 1e3
        match, dr2 = bestMatch(itau_vis, event.taus)
        
        if dr2<0.2*0.2:
            reco_tau = match

        # reco tau
        reco_tau_bis = None
        match = None
        dr2 = 1e3
        match, dr2 = bestMatch(itau_vis, event.taus_bis)
        
        if dr2<0.2*0.2:
            reco_tau_bis = match
    
        # identify the primary vertex
        event.the_pv = event.the_hn.vertex()

        # identify the secondary vertex
        event.the_sv = itau.vertex()

        # 2D transverse and 3D displacement, Pythagoras
        event.Lxy  = np.sqrt((event.the_pv.x() - event.the_sv.x())**2 + \
                             (event.the_pv.y() - event.the_sv.y())**2)

        event.Lxyz = np.sqrt((event.the_pv.x() - event.the_sv.x())**2 + \
                             (event.the_pv.y() - event.the_sv.y())**2 + \
                             (event.the_pv.z() - event.the_sv.z())**2)

        # reset before filling
        for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling    

        tofill['run'  ] = event.eventAuxiliary().run()
        tofill['lumi' ] = event.eventAuxiliary().luminosityBlock()
        tofill['event'] = event.eventAuxiliary().event()
     
        tofill['gen_tau_pt'   ] = itau.pt()     
        tofill['gen_tau_eta'  ] = itau.eta()    
        tofill['gen_tau_phi'  ] = itau.phi()    
        tofill['gen_tau_mass' ] = itau.mass()   
        tofill['gen_tau_pdgid'] = itau.pdgId()   
        tofill['gen_tau_q'    ] = itau.charge()   
        tofill['gen_tau_dm'   ] = itau.gen_dm  

        tofill['gen_tau_vis_pt'  ] = itau_vis.pt()     
        tofill['gen_tau_vis_eta' ] = itau_vis.eta()    
        tofill['gen_tau_vis_phi' ] = itau_vis.phi()    
        tofill['gen_tau_vis_mass'] = itau_vis.mass()   

        if reco_tau is not None:
            tofill['reco_tau_pt'   ] = reco_tau.pt()     
            tofill['reco_tau_eta'  ] = reco_tau.eta()  
            tofill['reco_tau_phi'  ] = reco_tau.phi()  
            tofill['reco_tau_mass' ] = reco_tau.mass() 
            tofill['reco_tau_pdgid'] = reco_tau.pdgId() 
            tofill['reco_tau_q'    ] = reco_tau.charge()   
            tofill['reco_tau_dxy'  ] = reco_tau.dxy()   
            tofill['reco_tau_dz'   ] = reco_tau.leadChargedHadrCand().dz() # forse 
            tofill['reco_tau_dm'   ] = reco_tau.decayMode()  # decay mode
#             tofill['reco_tau_dtj'  ] = reco_tau.pdgId()  # deeptau 
#             tofill['reco_tau_dte'  ] = reco_tau.pdgId()  # deeptau 
#             tofill['reco_tau_dtm'  ] = reco_tau.pdgId()  # deeptau 
            tofill['reco_tau_idj'  ] = 3 * reco_tau.tauID('byTightIsolationMVArun2v1DBnewDMwLT')  + \
                                       2 * reco_tau.tauID('byMediumIsolationMVArun2v1DBnewDMwLT') + \
                                       1 * reco_tau.tauID('byLooseIsolationMVArun2v1DBnewDMwLT')
            tofill['reco_tau_ide'  ] = 3 * reco_tau.tauID('againstElectronTightMVA6') + \
                                       2 * reco_tau.tauID('againstElectronMediumMVA6') + \
                                       1 * reco_tau.tauID('againstElectronLooseMVA6')  
            tofill['reco_tau_idm'  ] = 3 * reco_tau.tauID('againstMuonTight3') + \
                                       1 * reco_tau.tauID('againstMuonLoose3')  

        if reco_tau_bis is not None:
            tofill['reco_tau_bis_pt'   ] = reco_tau_bis.pt()     
            tofill['reco_tau_bis_eta'  ] = reco_tau_bis.eta()  
            tofill['reco_tau_bis_phi'  ] = reco_tau_bis.phi()  
            tofill['reco_tau_bis_mass' ] = reco_tau_bis.mass() 
            tofill['reco_tau_bis_pdgid'] = reco_tau_bis.pdgId() 
            tofill['reco_tau_bis_q'    ] = reco_tau_bis.charge()   
            tofill['reco_tau_bis_dxy'  ] = reco_tau_bis.dxy()   
            tofill['reco_tau_bis_dz'   ] = reco_tau_bis.leadChargedHadrCand().dz() # forse 
            tofill['reco_tau_bis_dm'   ] = reco_tau_bis.decayMode()  # decay mode
#             tofill['reco_tau_bis_dtj'  ] = reco_tau_bis.pdgId()  # deeptau 
#             tofill['reco_tau_bis_dte'  ] = reco_tau_bis.pdgId()  # deeptau 
#             tofill['reco_tau_bis_dtm'  ] = reco_tau_bis.pdgId()  # deeptau 
#             tofill['reco_tau_bis_idj'  ] = 3 * reco_tau_bis.tauID('byTightIsolationMVArun2v1DBnewDMwLT')  + \
#                                        2 * reco_tau_bis.tauID('byMediumIsolationMVArun2v1DBnewDMwLT') + \
#                                        1 * reco_tau_bis.tauID('byLooseIsolationMVArun2v1DBnewDMwLT')
#             tofill['reco_tau_bis_ide'  ] = 3 * reco_tau_bis.tauID('againstElectronTightMVA6') + \
#                                        2 * reco_tau_bis.tauID('againstElectronMediumMVA6') + \
#                                        1 * reco_tau_bis.tauID('againstElectronLooseMVA6')  
#             tofill['reco_tau_bis_idm'  ] = 3 * reco_tau_bis.tauID('againstMuonTight3') + \
#                                        1 * reco_tau_bis.tauID('againstMuonLoose3')  

        tofill['Lxy' ] = event.Lxy
        tofill['Lxyz'] = event.Lxyz
    
        ntuple.Fill(array('f',tofill.values()))

outfile.cd()
ntuple.Write()
outfile.Close()
