import ROOT
import pandas
import numpy as np
import root_numpy

import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
		    description="Convert MiniAOD to flat ntuples!")
parser.add_argument(
	"--file",
	choices=['ZTT','QCD'],
	required=True,
	help='Specify the sample you want to flatten')
args = parser.parse_args()
sample = args.file

# root_file = ROOT.TFile.Open('tau_gentau_tuple_{}.root'.format(sample), 'READ')
root_file = ROOT.TFile.Open('tau_gentau_tuple_taugun.root', 'READ')
tree = root_file.Get('tree')

dataframe = pandas.DataFrame( root_numpy.tree2array(    tree      = tree,
                                                        branches  = ['tau_gen_vis_signal_dR', 'tau_gen_vis_pt'],
                                                        selection = 'tau_gen_decaymode <= 10'))

QTL = [0.8, 0.85, 0.9, 0.95, 0.99]
pt_bins = [18., 20., 22., 24., 28., 32., 36., 40., 45., 50., 60., 80., 100., 200]
pt_bins_center = [0.5*pt_bins[jj] + 0.5*pt_bins[jj+1] for jj in range(len(pt_bins) - 1 )]

for ii, qq in enumerate(QTL):
	qtl = np.array([0.] * (len(pt_bins) - 1))
	for jj in range(len(pt_bins) - 1):
		selected = dataframe[dataframe.tau_gen_vis_pt >= pt_bins[jj]	]
		selected = selected [selected .tau_gen_vis_pt <  pt_bins[jj+1]	]
		if ii == 0: 
			print '[{}, {}]'.format(pt_bins[jj], pt_bins[jj+1]) , len(selected)
		qtl[jj]  = selected['tau_gen_vis_signal_dR'].quantile(qq)

	plt.plot(pt_bins_center, qtl, label = qq)

plt.xlabel('tau_gen_vis_pt [GeV]')
plt.ylabel('tau_gen_vis_signal_dR')
plt.xscale('log')

plt.grid()
plt.legend()
plt.savefig('quantile.pdf')
