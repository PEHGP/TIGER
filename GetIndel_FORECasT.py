#!/usr/bin/env python
import os,re,sys,collections,pandas,glob
from Bio.Seq import Seq
def get_pam(grna,ref):
	m=re.search(grna,ref)
	if m:
		start,end=m.span()
	else:
		ref=str(Seq(ref).reverse_complement())
		#print(grna)
		#print(ref)
		m=re.search(grna,ref)
		start,end=m.span()
	return ref,end
def get_indel_frequency(indelsummary):
	d=collections.defaultdict(int)
	for x in open(indelsummary):
		x=x.rstrip()
		l=x.split("\t")
		if l[0]=="-":
			continue
		indel=l[0].split("_")[0]
		if indel.startswith("I"):
			indel="+"+indel.split("I")[1]
		elif indel.startswith("D"):
			indel="-"+indel.split("D")[1]
		d[indel]+=int(l[2])
	total=sum(d.values())
	for i in d:
		d[i]=d[i]/float(total)*100
	return d
if __name__ == '__main__':
	flist=glob.glob("*_indel.xls")
	fd=collections.defaultdict(list)
	for f in flist:
		grna=f.split("_")[1]
		fd[grna].append(f)
	for x in open("200709-inDelphi-input.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		print(x)
		if len(fd[l[1]])==2:
			continue
		wt_prefix="wt_"+l[1]
		mut_prefix="mut_"+l[1]
		wt_grna=l[2].upper()
		mut_grna=l[3].upper()
		wt_ref=l[6].upper()
		mut_ref=l[7].upper()
		wt_ref,wt_pam=get_pam(wt_grna,wt_ref)
		mut_ref,mut_pam=get_pam(mut_grna,mut_ref)
		wt_ngg=wt_ref[wt_pam:wt_pam+3]
		mut_ngg=mut_ref[mut_pam:mut_pam+3]
		if (not wt_ngg.endswith("GG")) or (not mut_ngg.endswith("GG")):
			continue
		#print(mut_ref[mut_pam:mut_pam+3])
		#os.system("FORECasT.py %s %s %s"%(wt_ref,wt_pam,wt_prefix))
		#os.system("FORECasT.py %s %s %s"%(mut_ref,mut_pam,mut_prefix))
		wt_r=get_indel_frequency("%s_predictedindelsummary.txt"%wt_prefix)
		mut_r=get_indel_frequency("%s_predictedindelsummary.txt"%mut_prefix)
		#print(wt_r)
		wt_df=pandas.DataFrame.from_dict(wt_r,orient='index',columns=["Predicted frequency"])
		wt_df["indel_num"]=wt_df.index.astype(int)
		mut_df=pandas.DataFrame.from_dict(mut_r,orient='index',columns=["Predicted frequency"])
		mut_df["indel_num"]=mut_df.index.astype(int)
		wt_df.sort_values(by='indel_num',ascending=False).loc[:,["Predicted frequency"]].to_csv("%s_indel.xls"%wt_prefix,sep="\t",index_label="Indel length")
		mut_df.sort_values(by='indel_num',ascending=False).loc[:,["Predicted frequency"]].to_csv("%s_indel.xls"%mut_prefix,sep="\t",index_label="Indel length")
