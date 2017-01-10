from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram,cosh,append,savez,load,any,argsort,histogram2d,random,insert
import os
import numpy
import json
from optparse import OptionParser
from sys import stdout,argv
import ROOT as r
import pdb
#from helper_functions import *

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
#plt.style.use('atlas')
import matplotlib.mlab as mlab

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="Directory containing input files",type=str, default=".")
parser.add_option("--inputFile", help="Input file name",type=str, default="out")
parser.add_option("--submitDir", help="Directory containing output files",type=str, default="output")
parser.add_option("--plotDir", help="Directory containing plots",type=str, default="plots")
parser.add_option("--numEvents", help="How many events to include (set to -1 for all events)",type=int, default=-1)
parser.add_option("--topMassList", help="Truth top masses",type=str, default="allTopMasses.json")
parser.add_option("-r","--root", help="force reading in root files",action="store_true", default=False)
parser.add_option("-s","--smear", help="force resmearing jets",action="store_true", default=False)
parser.add_option("-l","--learn", help="force relearning jets",action="store_true", default=False)

# Root configuration
parser.add_option("--j1_pt", help="jet 1 pT branch name",type=str, default="j1_pt")
parser.add_option("--j1_eta", help="jet 1 eta branch name",type=str, default="j1_eta")
parser.add_option("--j1_phi", help="jet 1 phi branch name",type=str, default="j1_phi")
parser.add_option("--j1_m", help="jet 1 m branch name",type=str, default="j1_m")
parser.add_option("--j2_pt", help="jet 2 pT branch name",type=str, default="j2_pt")
parser.add_option("--j2_eta", help="jet 2 eta branch name",type=str, default="j2_eta")
parser.add_option("--j2_phi", help="jet 2 phi branch name",type=str, default="j2_phi")
parser.add_option("--j2_m", help="jet 2 m branch name",type=str, default="j2_m")
parser.add_option("--b_pt", help="jet 1 pT branch name",type=str, default="b_pt")
parser.add_option("--b_eta", help="jet 1 eta branch name",type=str, default="b_eta")
parser.add_option("--b_phi", help="jet 1 phi branch name",type=str, default="b_phi")
parser.add_option("--b_m", help="jet 1 m branch name",type=str, default="b_m")

(options, args) = parser.parse_args()

if not os.path.exists(options.inputDir): raise OSError(options.inputDir +' does not exist. This is where the input Root files go.')
if not os.path.exists(options.submitDir):
  print '== Making folder '+options.submitDir+' =='
  os.makedirs(options.submitDir)
if not os.path.exists(options.plotDir):
  print '== Making folder '+options.plotDir+' =='
  os.makedirs(options.plotDir)
if not os.path.exists(options.topMassList): raise OSError(options.topMassList +' does not exist. This is where the list of top masses go.')

def read(topMass):
  import glob
  global cutflow

  if len(options.inputFile)>0:
    filename = options.inputFile+topMass+".root"
    filenames = glob.glob(options.inputDir+'/'+filename)
    if len(filenames) == 0: raise OSError('Can\'t find file '+options.inputDir+'/'+options.inputFile) 
  else:
    filenames = glob.glob(options.inputDir+'/*.root')
    if len(filenames) == 0: raise OSError('Can\'t find files in '+options.inputDir) 
  tree = r.TChain('oTree')
  for filename in filenames:
    statinfo = os.stat(filename)
    if statinfo.st_size < 10000: continue #sometimes batch jobs fail
    print '== Reading in '+filename+' =='

    tree.Add(filename) 

  # make sure the branches are compatible between the two
  branches = set(i.GetName() for i in tree.GetListOfBranches())

  # required:
  if options.j1_pt not in branches: raise RuntimeError(options.j1_pt+' branch does not exist. This is the branch containing jet 1 pTs.')
  else: print '== \''+options.j1_pt+'\' branch is being read as jet 1 pTs =='
  if options.j1_eta not in branches: raise RuntimeError(options.j1_eta+' branch does not exist. This is the branch containing jet 1 etas.')
  else: print '== \''+options.j1_eta+'\' branch is being read as jet 1 etas =='
  if options.j1_phi not in branches: raise RuntimeError(options.j1_phi+' branch does not exist. This is the branch containing jet 1 phis.')
  else: print '== \''+options.j1_phi+'\' branch is being read as jet 1 phis =='
  if options.j1_m not in branches: raise RuntimeError(options.j1_m+' branch does not exist. This is the branch containing jet 1 ms.')
  else: print '== \''+options.j1_m+'\' branch is being read as jet 1 ms =='
  if options.j2_pt not in branches: raise RuntimeError(options.j2_pt+' branch does not exist. This is the branch containing jet 2 pTs.')
  else: print '== \''+options.j2_pt+'\' branch is being read as jet 2 pTs =='
  if options.j2_eta not in branches: raise RuntimeError(options.j2_eta+' branch does not exist. This is the branch containing jet 2 etas.')
  else: print '== \''+options.j2_eta+'\' branch is being read as jet 2 etas =='
  if options.j2_phi not in branches: raise RuntimeError(options.j2_phi+' branch does not exist. This is the branch containing jet 2 phis.')
  else: print '== \''+options.j2_phi+'\' branch is being read as jet 2 phis =='
  if options.j2_m not in branches: raise RuntimeError(options.j2_m+' branch does not exist. This is the branch containing jet 2 ms.')
  else: print '== \''+options.j2_m+'\' branch is being read as jet 2 ms =='
  if options.b_pt not in branches: raise RuntimeError(options.b_pt+' branch does not exist. This is the branch containing b jet pTs.')
  else: print '== \''+options.b_pt+'\' branch is being read as b jet pTs =='
  if options.b_eta not in branches: raise RuntimeError(options.b_eta+' branch does not exist. This is the branch containing b jet etas.')
  else: print '== \''+options.b_eta+'\' branch is being read as b jet etas =='
  if options.b_phi not in branches: raise RuntimeError(options.b_phi+' branch does not exist. This is the branch containing b jet phis.')
  else: print '== \''+options.b_phi+'\' branch is being read as b jet phis =='
  if options.b_m not in branches: raise RuntimeError(options.b_m+' branch does not exist. This is the branch containing b jet ms.')
  else: print '== \''+options.b_m+'\' branch is being read as b jet ms =='

  nentries = tree.GetEntries()

  t_jet1s = []
  t_jet2s = []
  t_bjets = []
  t_Ws = []
  t_ts = []

  for jentry in xrange(nentries):
    if jentry>options.numEvents and options.numEvents>0: continue
    tree.GetEntry(jentry)

    if not jentry%1000:
      stdout.write('== \r%d events read ==\n'%jentry)
      stdout.flush() 

    #Ntruthphotons = getattr(tree,options.Ntruthphotons)
    #if not Ntruthphotons == 2: continue #only events with 2 truth photons

    jet1 = r.TLorentzVector()
    jet1.SetPtEtaPhiM(getattr(tree,options.j1_pt),getattr(tree,options.j1_eta),getattr(tree,options.j1_phi),getattr(tree,options.j1_m))
    t_jet1s.append(jet1)

    jet2 = r.TLorentzVector()
    jet2.SetPtEtaPhiM(getattr(tree,options.j2_pt),getattr(tree,options.j2_eta),getattr(tree,options.j2_phi),getattr(tree,options.j2_m))
    t_jet2s.append(jet2)

    bjet = r.TLorentzVector()
    bjet.SetPtEtaPhiM(getattr(tree,options.b_pt),getattr(tree,options.b_eta),getattr(tree,options.b_phi),getattr(tree,options.b_m))
    t_bjets.append(bjet)

    W = jet1+jet2
    t_Ws.append(W)
    t_ts.append(W+bjet)

  return t_jet1s,t_jet2s,t_bjets,t_Ws,t_ts

def smear(jets):
  smeared_jets = []
  smears = random.normal(0,7,len(jets)) #smear jets by 7 GeV
  for j,s in zip(jets,smears):
    sj = r.TLorentzVector()
    sj.SetPtEtaPhiM(j.Pt()+s,j.Eta(),j.Phi(),j.M())
    smeared_jets.append(sj)
  return array(smeared_jets)

def learn(jets,tjets,label='',name=''):
  X = array([j.Pt() for j in jets])
  Y = array([tj.Pt() for tj in tjets])
  model = Sequential()
  model.add(Dense(5, input_dim=1, init='uniform', activation='sigmoid'))
  #model.add(Dense(8, init='uniform', activation='relu'))
  model.add(Dense(1, init='uniform')) # regression
  model.compile(loss='mean_squared_error', optimizer='adam', metrics=[])
  model.fit(X, Y, nb_epoch=30, batch_size=10)
  L = model.predict(X)
  L = array([L[i][0] for i in range(len(L))]) #reshape

  learned_jets = []
  for lpt,j in zip(L,jets):
    lj = r.TLorentzVector()
    lj.SetPtEtaPhiM(lpt,j.Eta(),j.Phi(),j.M())
    learned_jets.append(lj)

  binwidth = 2
  H, xedges, yedges = histogram2d(Y,X,bins=[numpy.arange(0, 100, binwidth),numpy.arange(0, 100, binwidth)])
  plt.pcolormesh(xedges,yedges,H.T,cmap=plt.get_cmap('Blues'))
  sampleX = numpy.arange(0,100,1)
  sampleY = model.predict(sampleX)
  plt.plot(sampleY,sampleX,ls='-',color='r')
  plt.xlabel(label+r' $p_T^{true}$')
  plt.ylabel(label+r' $p_T^{reco}$')
  plt.xlim(0,100)
  plt.ylim(0,100)
  plotname = name+'_pt_hist2d'
  plt.savefig(options.plotDir+'/'+plotname+'.png')
  plt.close()
  
  return array(learned_jets)

def plot(objects,label='',name=''):
  colors = ['b','r','black','g','orange','purple','yellow']

  binwidth = 2
  for topMass,color in zip(objects.keys(),colors):
    data = [o.Pt() for o in objects[topMass]]
    n,bins = numpy.histogram(data,normed=True,bins=numpy.arange(0, 100, binwidth))
    n = insert(n,0,0)
    n = insert(n,len(n),0)
    bins = insert(bins,len(bins),bins[-1]+binwidth)
    n = n/sum(n)
    plt.plot(bins,n,drawstyle='steps',ls='-',color=color,label='$M_t=$'+topMass)
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,100)
  plt.xlabel(label+r' $p_T$')
  plt.ylabel('a.u.')
  plt.legend(loc='upper right',frameon=False)
  plotname = name+'_pt'
  plt.savefig(options.plotDir+'/'+plotname+'.png')
  plt.close()

  for topMass,color in zip(objects.keys(),colors):
    data = [o.Eta() for o in objects[topMass]]
    binwidth = 0.1
    n,bins = numpy.histogram(data,normed=True,bins=numpy.arange(-3, 3, binwidth))
    n = insert(n,0,0)
    n = insert(n,len(n),0)
    bins = insert(bins,len(bins),bins[-1]+binwidth)
    n = n/sum(n)
    plt.plot(bins,n,drawstyle='steps',ls='-',color=color,label='$M_t=$'+topMass)
  plt.ylim(0,max(n)*1.1)
  plt.xlim(-3,3)
  plt.xlabel(label+r' $\eta$')
  plt.ylabel('a.u.')
  plt.legend(loc='upper right',frameon=False)
  plotname = name+'_eta'
  plt.savefig(options.plotDir+'/'+plotname+'.png')
  plt.close()

  for topMass,color in zip(objects.keys(),colors):
    data = [o.Phi() for o in objects[topMass]]
    binwidth = 0.1
    n,bins = numpy.histogram(data,normed=True,bins=numpy.arange(-3.2, 3.2, binwidth))
    n = insert(n,0,0)
    n = insert(n,len(n),0)
    bins = insert(bins,len(bins),bins[-1]+binwidth)
    n = n/sum(n)
    plt.plot(bins,n,drawstyle='steps',ls='-',color=color,label='$M_t=$'+topMass)
  plt.ylim(0,max(n)*1.1)
  plt.xlim(-3.2,3.2)
  plt.xlabel(label+r' $\phi$')
  plt.ylabel('a.u.')
  plt.legend(loc='upper right',frameon=False)
  plotname = name+'_phi'
  plt.savefig(options.plotDir+'/'+plotname+'.png')
  plt.close()

  for topMass,color in zip(objects.keys(),colors):
    data = [o.M() for o in objects[topMass]]
    print mean(data),std(data),name,topMass
    binwidth = 5 
    n,bins = numpy.histogram(data,normed=True,bins=numpy.arange(0, 250, binwidth))
    n = insert(n,0,0)
    n = insert(n,len(n),0)
    bins = insert(bins,len(bins),bins[-1]+binwidth)
    n = n/sum(n)
    plt.plot(bins,n,drawstyle='steps',ls='-',color=color,label='$M_t=$'+topMass)
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,250)
  plt.xlabel(label+r' $M$')
  plt.ylabel('a.u.')
  plt.legend(loc='upper right',frameon=False)
  plotname = name+'_m'
  plt.savefig(options.plotDir+'/'+plotname+'.png')
  plt.close()

with open(options.topMassList) as data_file:    
  topMasses = json.load(data_file)
  topMasses = [str(tm) for tm in topMasses]

t_jet1s = {m:[] for m in topMasses}
t_jet2s = {m:[] for m in topMasses}
t_bjets = {m:[] for m in topMasses}
t_Ws = {m:[] for m in topMasses}
t_ts = {m:[] for m in topMasses}
jet1s = {m:[] for m in topMasses}
jet2s = {m:[] for m in topMasses}
bjets = {m:[] for m in topMasses}
Ws = {m:[] for m in topMasses}
ts = {m:[] for m in topMasses}
l_jet1s = {m:[] for m in topMasses}
l_jet2s = {m:[] for m in topMasses}
l_bjets = {m:[] for m in topMasses}
l_Ws = {m:[] for m in topMasses}
l_ts = {m:[] for m in topMasses}

for topMass in topMasses:
  doRoot = options.root
  if not doRoot:
    try:
      t_jet1s[topMass] = load(options.submitDir+'/'+'t_jet1s'+'_'+topMass+'.npy')
      t_jet2s[topMass] = load(options.submitDir+'/'+'t_jet2s'+'_'+topMass+'.npy')
      t_bjets[topMass] = load(options.submitDir+'/'+'t_bjets'+'_'+topMass+'.npy')
      t_Ws[topMass] = load(options.submitDir+'/'+'t_Ws'+'_'+topMass+'.npy')
      t_ts[topMass] = load(options.submitDir+'/'+'t_ts'+'_'+topMass+'.npy')
    except IOError:
      print 'Numpy files for truth data don\'t exist.'
      doRoot = True
  if doRoot:
    'Attempting to read in Root files.'
    t_jet1s[topMass],t_jet2s[topMass],t_bjets[topMass],t_Ws[topMass],t_ts[topMass] = read(topMass)

    save(options.submitDir+'/'+'t_jet1s'+'_'+topMass+'.npy',t_jet1s[topMass])
    save(options.submitDir+'/'+'t_jet2s'+'_'+topMass+'.npy',t_jet2s[topMass])
    save(options.submitDir+'/'+'t_bjets'+'_'+topMass+'.npy',t_bjets[topMass])
    save(options.submitDir+'/'+'t_Ws'+'_'+topMass+'.npy',t_Ws[topMass])
    save(options.submitDir+'/'+'t_ts'+'_'+topMass+'.npy',t_ts[topMass])

  doSmear = options.smear
  if not doSmear:
    try:
      jet1s[topMass] = load(options.submitDir+'/'+'jet1s'+'_'+topMass+'.npy')
      jet2s[topMass] = load(options.submitDir+'/'+'jet2s'+'_'+topMass+'.npy')
      bjets[topMass] = load(options.submitDir+'/'+'bjets'+'_'+topMass+'.npy')
      Ws[topMass] = load(options.submitDir+'/'+'Ws'+'_'+topMass+'.npy')
      ts[topMass] = load(options.submitDir+'/'+'ts'+'_'+topMass+'.npy')
    except IOError:
      print 'Numpy files for smeared data don\'t exist.'
      doSmear = True
  if doSmear:
    print 'Smearing.'
    jet1s[topMass] = smear(t_jet1s[topMass])
    jet2s[topMass] = smear(t_jet2s[topMass])
    bjets[topMass] = smear(t_bjets[topMass])
    Ws[topMass] = jet1s[topMass]+jet2s[topMass]
    ts[topMass] = Ws[topMass]+bjets[topMass]

    save(options.submitDir+'/'+'jet1s'+'_'+topMass+'.npy',jet1s[topMass])
    save(options.submitDir+'/'+'jet2s'+'_'+topMass+'.npy',jet2s[topMass])
    save(options.submitDir+'/'+'bjets'+'_'+topMass+'.npy',bjets[topMass])
    save(options.submitDir+'/'+'Ws'+'_'+topMass+'.npy',Ws[topMass])
    save(options.submitDir+'/'+'ts'+'_'+topMass+'.npy',ts[topMass])

  doLearn = options.learn
  if not doLearn:
    try:
      l_jet1s[topMass] = load(options.submitDir+'/'+'l_jet1s'+'_'+topMass+'.npy')
      l_jet2s[topMass] = load(options.submitDir+'/'+'l_jet2s'+'_'+topMass+'.npy')
      l_bjets[topMass] = load(options.submitDir+'/'+'l_bjets'+'_'+topMass+'.npy')
      l_Ws[topMass] = load(options.submitDir+'/'+'l_Ws'+'_'+topMass+'.npy')
      l_ts[topMass] = load(options.submitDir+'/'+'l_ts'+'_'+topMass+'.npy')
    except IOError:
      print 'Numpy files for learned data don\'t exist.'
      doLearn = True
  if doLearn:
    print 'Learning.'
    from keras.models import Sequential
    from keras.layers import Dense
    from sklearn.preprocessing import StandardScaler
    l_jet1s[topMass] = learn(jet1s[topMass],t_jet1s[topMass],label='Jet 1',name='jet1'+'_'+topMass)
    l_jet2s[topMass] = learn(jet2s[topMass],t_jet2s[topMass],label='Jet 2',name='jet2'+'_'+topMass)
    l_bjets[topMass] = learn(bjets[topMass],t_bjets[topMass],label='b-Jet',name='bjet'+'_'+topMass)
    l_Ws[topMass] = l_jet1s[topMass]+l_jet2s[topMass]
    l_ts[topMass] = l_Ws[topMass]+l_bjets[topMass]

    save(options.submitDir+'/'+'l_jet1s'+'_'+topMass+'.npy',l_jet1s[topMass])
    save(options.submitDir+'/'+'l_jet2s'+'_'+topMass+'.npy',l_jet2s[topMass])
    save(options.submitDir+'/'+'l_bjets'+'_'+topMass+'.npy',l_bjets[topMass])
    save(options.submitDir+'/'+'l_Ws'+'_'+topMass+'.npy',l_Ws[topMass])
    save(options.submitDir+'/'+'l_ts'+'_'+topMass+'.npy',l_ts[topMass])

plot(t_jet1s,label='Jet 1',name='tjet1')
plot(t_jet2s,label='Jet 2',name='tjet2')
plot(t_bjets,label='b-Jet',name='tbjet')
plot(t_Ws,label='W',name='tW')
plot(t_ts,label='Top',name='tt')

plot(jet1s,label='Jet 1',name='jet1')
plot(jet2s,label='Jet 2',name='jet2')
plot(bjets,label='b-Jet',name='bjet')
plot(Ws,label='W',name='W')
plot(ts,label='Top',name='t')

plot(l_jet1s,label='Jet 1',name='ljet1')
plot(l_jet2s,label='Jet 2',name='ljet2')
plot(l_bjets,label='b-Jet',name='lbjet')
plot(l_Ws,label='W',name='lW')
plot(l_ts,label='Top',name='lt')
