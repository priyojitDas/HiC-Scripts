import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from pylab import rcParams
plt.rcParams.update({'font.size': 25})
rcParams['figure.figsize'] = 20, 10
import sys

def input_nd_process_hic(loc):
	df = pd.read_csv(loc,comment='#',sep='\t')
	df = df.set_index(list(df)[0])
	#df = df.dropna(axis='columns',how='all')
	#df = df.dropna(axis='index',how='all')	
	df = df.fillna(0)
	indx = list(df)
	return df.values, indx

c1 = sys.argv[1]
c2 = sys.argv[2]
c3 = sys.argv[3]

cell_list = [c1]
reps = [c2]
bases = [c3]

zval_arr = []

chrm = [_ for _ in range(1,23)] + ['X']

pos = []
strength = []
istrength = []

for ch in chrm:
	print ch
	comp_strn = []
	indv_strn = []
	for cindx,cell in enumerate(cell_list):
		zdf,indx = input_nd_process_hic(bases[cindx]+str(ch)+".zScore.matrix.gz")
		dcomp = np.nan_to_num(np.array(pd.read_csv(bases[cindx]+str(ch)+".zScore.eigen1.bedGraph",header=None,sep='\t',comment='t')[3].astype('float')))
		
		dcomp_zero = np.where(dcomp == 0)[0]

		dcomp_filtered = np.delete(dcomp,dcomp_zero)

		zdf_filtered = np.delete(zdf,dcomp_zero,axis=0)
		zdf_filtered = np.delete(zdf_filtered,dcomp_zero,axis=1)

		#sns.heatmap(np.log2(zdf_filtered+1),cmap='RdBu',center=0)
		#plt.show()

		#exit()

		dcomp_filtered_sorted = np.argsort(dcomp_filtered)

		zdf_filtered_sorted = zdf_filtered[dcomp_filtered_sorted,:]

		zdf_filtered_sorted = zdf_filtered_sorted[:,dcomp_filtered_sorted]

		dcomp_filtered_sorted_binary = dcomp_filtered[dcomp_filtered_sorted]

		dcomp_filtered_sorted_binary[dcomp_filtered[dcomp_filtered_sorted] > 0] = 1
		dcomp_filtered_sorted_binary[dcomp_filtered[dcomp_filtered_sorted] <= 0] = 3

		#plt.plot(dcomp_filtered_sorted_binary)
		#plt.show()

		#sns.heatmap(np.dot(dcomp_filtered_sorted_binary.reshape((dcomp_filtered_sorted.shape[0],1)),dcomp_filtered_sorted_binary.reshape((1,dcomp_filtered_sorted.shape[0]))),cmap='seismic')
		#plt.show()
		
		comp_dot = np.dot(dcomp_filtered_sorted_binary.reshape((dcomp_filtered_sorted.shape[0],1)),dcomp_filtered_sorted_binary.reshape((1,dcomp_filtered_sorted.shape[0])))

		comp_dot = np.array(comp_dot,dtype='int32')

		#sns.heatmap(zdf_filtered_sorted,cmap='seismic')
		#plt.show()

		zdf_filtered_sorted_cg = []
		dcomp_filtered_sorted_cg = []

		for i in xrange(0,zdf_filtered_sorted.shape[0],2):
			tmpLst1 = []
			tmpLst2 = []
			for j in xrange(0,zdf_filtered_sorted.shape[0],2):
				tmpLst1.append(np.mean(zdf_filtered_sorted[i:i+2,j:j+2]))
				tmpLst2.append(np.bincount(comp_dot[i:i+2,j:j+2].flatten()).argmax())
			zdf_filtered_sorted_cg.append(tmpLst1)
			dcomp_filtered_sorted_cg.append(tmpLst2)

		zdf_filtered_sorted_cg = np.array(zdf_filtered_sorted_cg)

		dcomp_filtered_sorted_cg = np.array(dcomp_filtered_sorted_cg)

		zdf_filtered_sorted_cg = gaussian_filter(zdf_filtered_sorted_cg,sigma=1)

		uniq_dcomp = list(set(dcomp_filtered_sorted_cg.flatten()))

		#print ch,cell,uniq_dcomp

		#print dcomp_filtered_sorted_cg

		#ns.heatmap(dcomp_filtered_sorted_cg,cmap='seismic',center=0)
		#plt.show()
	
		comp_BB = zdf_filtered_sorted_cg[np.triu_indices(zdf_filtered_sorted_cg.shape[0])]
		comp_BB = comp_BB[dcomp_filtered_sorted_cg[np.triu_indices(dcomp_filtered_sorted_cg.shape[0])] == uniq_dcomp[0]]
		#comp_BB = zdf_filtered_sorted_cg[dcomp_filtered_sorted_cg == uniq_dcomp[0]]
		comp_tBB = np.median(comp_BB[comp_BB > np.percentile(comp_BB,80)])

		comp_BA = zdf_filtered_sorted_cg[np.triu_indices(zdf_filtered_sorted_cg.shape[0])]
		comp_BA = comp_BA[dcomp_filtered_sorted_cg[np.triu_indices(dcomp_filtered_sorted_cg.shape[0])] == uniq_dcomp[1]]
		#comp_tBA = np.median(comp_BA[comp_BA < np.percentile(comp_BA,20)])
		comp_tBA = np.median(comp_BA[comp_BA > np.percentile(comp_BA,80)])

		comp_AA = zdf_filtered_sorted_cg[np.triu_indices(zdf_filtered_sorted_cg.shape[0])]
                comp_AA = comp_AA[dcomp_filtered_sorted_cg[np.triu_indices(dcomp_filtered_sorted_cg.shape[0])] == uniq_dcomp[2]]
		#comp_AA = zdf_filtered_sorted_cg[dcomp_filtered_sorted_cg == uniq_dcomp[2]]
		comp_tAA = np.median(comp_AA[comp_AA > np.percentile(comp_AA,80)])

		#print cell,(comp_tAA + comp_tBB) / (2 * comp_tBA)

		#comp_strn.append((comp_tAA + comp_tBB) / (2 * comp_tBA))
		comp_strn.append((comp_tAA + comp_tBB) - (2 * comp_tBA))
		indv_strn.append([comp_tBB-comp_tBA,comp_tAA-comp_tBA])

		#sns.heatmap(zdf_filtered_sorted_cg,cmap='seismic',center=0)
		#plt.show()

		indx_filtered = np.delete(indx,dcomp_zero)[:zdf_filtered_sorted_cg.shape[0]]

		#sns.heatmap(zdf_filtered_sorted_cg,cmap='seismic',center=0)
		#plt.show()

		#zmean = zdf_filtered_sorted_cg.mean(axis=1)
		#zval_arr.append(list(zmean))
		
		#plt.bar([cindx*2,cindx*2+1],[comp_tBB-comp_tBA,comp_tAA-comp_tBA],facecolor=color[cindx],label=cell_list[cindx],width = 0.4,edgecolor='black')
	pos.append(list(np.argsort(comp_strn)))
	strength.append(list(comp_strn))
	istrength.append(indv_strn)
	#plt.plot([_*2+0.5 for _ in xrange(len(cell_list))],comp_strn,linestyle='--',color='black',linewidth=3,marker='X')
	#plt.xticks([_*2+0.5 for _ in xrange(len(cell_list))],scell_list)
	#plt.ylabel('Rank')
	#plt.savefig('Saddle1Cohort/chr'+str(ch)+'.png',dpi=300)
	#plt.clf()
	#plt.close()

pos = np.array(pos)
strength = np.array(strength)
istrength = np.array(istrength)
print strength
print istrength

for cindx,cell in enumerate(cell_list):
	fig, ax1 = plt.subplots()
	for idx,ch in enumerate(chrm):
		ax1.bar([idx*2+1-0.2,idx*2+1+0.2],istrength[idx,cindx,:],color=['blue','red'],width=0.2)#facecolor=color[cindx],width = 0.2,edgecolor=color[cindx])
	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	ax2.plot([_*2+1 for _ in xrange(len(chrm))],strength[:,cindx],linestyle='--',color='black',linewidth=3,marker='X',label=cell_list[0])
	plt.xticks([_*2+1 for _ in xrange(len(chrm))],chrm)
	ax1.set_ylabel('Compartment Strength')
	ax1.set_ylim([0,6])
	ax2.set_ylim([1,8])
	ax2.set_ylabel('Average Compartment Strength')
	#ax2.legend(loc=1)
	plt.savefig(cell_list[cindx]+'_'+reps[cindx]+'_S1v2.png',dpi=300,bbox_inches='tight', pad_inches=0)
	plt.clf()
	plt.close()