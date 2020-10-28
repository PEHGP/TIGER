#!/usr/bin/env python
import matplotlib,re,glob,t7,os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import random,pandas,sys,itertools
from Bio import SeqIO
from matplotlib.colors import ListedColormap,BoundaryNorm
import stat_many
def get_consurf(consurf_file,wt_align_seq):
	cl=[]
	seq=""
	for x in open(consurf_file):
		x=x.rstrip()
		l=x.split("\t")
		cl.append(int(l[1]))
		seq+=l[0]
	#print seq
	#print wt_align_seq.replace("-","")
	m=re.search(wt_align_seq.replace("-",""),seq)
	start,end=m.span()
	cl=cl[start:end]
	consurfl=[]
	t=0
	for s in wt_align_seq:
		if s!="-":
			consurfl.append(cl[t])
			t+=1
		else:
			consurfl.append(0)
	return consurfl
def get_data(prefix,consurf_file,df_summary):
	top=40
	df=pandas.read_csv("%s_alleles_deal.xls"%prefix,sep='\t',header=0,index_col=0)
	df_group=pandas.concat([df,df_summary],axis=1,join='inner')
	#print(df_group)
	df_func=df_group.groupby("AA_seq")["p_function"].max().to_frame(name="p_function")
	df_func2=df_group.groupby("AA_seq")["p_function"].min().to_frame(name="p_function")
	for index,row in df_func.iterrows():
		if(row["p_function"]!=df_func2.loc[index,"p_function"]):
			print("seq function is not consistent")
			sys.exit()
	df_sum=df_group.groupby("AA_seq")["#Reads"].sum().to_frame(name="#Reads")
	df_sum["%Reads"]=df_sum["#Reads"]/df_sum["#Reads"].sum()*100
	df_sum=pandas.concat([df_sum,df_func],axis=1,join='inner')
	df_sum.sort_values(by=["#Reads"],ascending=False,inplace=True)
	print(df_sum)
	wt_seq=df_summary.loc["wt","AA_seq"]
	#dseq=SeqIO.index("%s_align.faa"%prefix,"fasta")
	consurfl=get_consurf(consurf_file,wt_seq)
	wt_seq=list(wt_seq)
	d={}
	t=0
	for seq in df_sum.head(top).index:
		#print(seq)
		d[seq]=list(seq)
		t+=1
	d["wt_cds"]=wt_seq
	df=pandas.DataFrame.from_dict(d,orient='index')
	#print (df!="-").any(axis=0)
	df=df.loc[:, (df != "-").any(axis=0)]
	df2=df.loc[["wt_cds"]+list(df_sum.head(top).index),:]
	print(df2)
	conser_score=[]
	indexd={}
	for name in df2.index:
		if name=="wt_cds":
			indexd[name]="wt"
		else:
			indexd[name]="%.2f%%(%s)%s"%(df_sum.loc[name,"%Reads"],df_sum.loc[name,"#Reads"],df_sum.loc[name,"p_function"])
		ri=[]
		for j in df2.columns:
			if df2.loc[name,j]=="-":
				if df2.loc["wt_cds",j]=="-":
					ri.append(-1)
				else:
					ri.append(0)
			else:
				if df2.loc["wt_cds",j]=="-":
					ri.append(-1)
				else:
					ri.append(int(consurfl[j]))
		conser_score.append(ri)
	df2.rename(index=indexd,inplace=True)
	np.savetxt("%s_heatmap_score.txt"%prefix,conser_score,delimiter="\t",fmt="%s")
	df2.to_csv("%s_heatmap_aa.xls"%prefix,sep="\t",index_label="reads")
	return df2,np.array(conser_score)
def get_heatmap(df,conser_score,prefix):
	print(df.shape)
	print(conser_score.shape)
	x,y=df.shape
	#fig = plt.figure(dpi=300)
	fig=plt.figure(figsize=(y*0.4,(x+1)*0.6),dpi=300)
	gs=matplotlib.gridspec.GridSpec(1,2,width_ratios=[1,30]	)
	#ax = fig.add_subplot(2,1,1)
	ax = plt.subplot(gs[1])
	#print(conser_score)
	#ax.imshow(conser_score,cmap='PiYG_r',aspect='equal')
	#masked_conser = np.ma.masked_where(conser_score ==0,conser_score)
	#cmap = plt.get_cmap('PiYG_r')
	#cmap.set_bad(color='grey')
	col=["#A2A2A2","#A2A2A2","#0A7D82","#4BAFBE","#A5DCE6","#D7F0F0","#FFFFFF","#FAEBF5","#FAC8DC","#F07DAA","#A0285F"]
	cold={}
	for t,i in enumerate(range(-1,10)):
		#print(t)
		cold[i]=col[t]
	print(col)
	col=[]
	sl=sorted(list(set(itertools.chain.from_iterable(conser_score))))
	print(sl)
	for i in sl:
		col.append(cold[i])
	print(col)
	#sys.exit()
	cmap=ListedColormap(col)
	print(conser_score)
	norm=BoundaryNorm(sl+[sl[-1]+1],cmap.N)
	h=ax.imshow(conser_score,cmap=cmap,aspect='equal',norm=norm)
	ax.set_yticks(np.arange(x))
	ax.set_yticklabels(list(df.index))
	ax.tick_params('y',left=False,right=True,labelleft=False,labelright=True)
	for edge, spine in ax.spines.items():
		spine.set_visible(False)
	ax.set_xticks(np.arange(df.shape[1])-.5, minor=True)
	ax.set_yticks(np.arange(df.shape[0]+1)-.5, minor=True)
	ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
	ax.tick_params(which="minor", bottom=False, left=False)
	for i in range(x):
		for j in range(y):
			#print i,j
			if df.iloc[i,j]!=df.loc["wt",df.columns[j]] and df.iloc[i,j]!="-":
				ax.text(j,i,df.iloc[i,j],ha="center",va="center",color="red")
			else:
				if df.iloc[i,j]=="-" and df.iloc[0,j]=="-":
					ax.text(j,i,"",ha="center",va="center")
				else:
					ax.text(j,i,df.iloc[i,j],ha="center",va="center")
			#if df.iloc[i,j]=="-":
			#	highlight_cell(j,i, color="grey",linewidth=3)
	ax.set_xticks([])
	#ax.set_yticks([])
	#plt.legend()
	#plt.colorbar(h,spacing='uniform')
	#plt.tight_layout()
	#ax2 = fig.add_subplot(2,1,2)
	ax2=plt.subplot(gs[0])
	col_leg=[]
	num_leg=[]
	for i in sl:
		if i==0 or i==-1:
			continue
		num_leg.append(i)
		col_leg.append(cold[i])
	print(col_leg)
	#sys.exit()
	cmap=ListedColormap(col_leg)
	print(np.array(num_leg))
	norm=BoundaryNorm(num_leg+[num_leg[-1]+1],cmap.N)
	ax2.imshow(np.array([num_leg]).T,cmap=cmap,aspect='equal',norm=norm)
	for i,n in enumerate(num_leg):
		ax2.text(0,i,str(n),ha="center",va="center")
	ax2.set_yticks([0,(len(num_leg)+1)/2-1,len(num_leg)-1])
	ax2.set_yticklabels(["Variable","Average","Conserved"])
	ax2.set_xticks([])
	#plt.subplots_adjust(wspace=0,hspace=0)
	fig.tight_layout()
	fig.savefig("%s_heatmap.png"%prefix,format='png')
	fig.clf()
	plt.close(fig)
if __name__ == '__main__':
	aa_file="mut_cds.fa"
	extend=20
	for x in open("input_051_3j.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		sl=os.popen("find -L . -name \"CRISPResso_on_%s\""%(l[0])).readlines()
		target_aa=l[-1].upper()
		wt_cds_seq=l[6].upper()
		for p in sl:
			p=p.rstrip()
			print(p)
			name=p.split("/")[1].split("-")[0]+"_"+l[0]
			sample=name.split("_")[1]
			if not os.path.exists("%s_alleles_deal.xls"%name):
				continue
			if os.path.exists("%s_heatmap.png"%name):
				continue
			#if name!="Rep2_56-051":
			#	continue
			grna_seq=t7.get_grna(l[5].upper(),sample)
			cut_point,orientation=t7.my_cut(grna_seq,wt_cds_seq)
			grna_prot=t7.get_grna_prot(wt_cds_seq,grna_seq,extend,cut_point)
			df_summary=stat_many.get_region(name+"_align.faa",grna_prot)
			consurf_file=t7.get_consurf_file(target_aa,aa_file)
			df,conser_score=get_data(name,consurf_file,df_summary)
			get_heatmap(df,conser_score,name)
			#break
		#break
