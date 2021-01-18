#!/usr/bin/env python
import re,pandas,glob,collections,sys,t7
from Bio import SeqIO
def get_region(align_fa,grna_prot):
	wt_name="wt_cds"
	seqd=SeqIO.index(align_fa,"fasta")
	wt_seq=str(seqd[wt_name].seq).upper()
	wt_seq_deal=wt_seq.replace("-","").upper()
	m=re.search(grna_prot,wt_seq_deal)
	start,end=m.span()
	pseq=[]
	for n,s in enumerate(wt_seq):
		if s!="-":
			pseq.append(n)
	real_start=pseq[start]
	if end==len(pseq):
		real_end=pseq[end-1]+1
	else:
		real_end=pseq[end]
	rd={}
	for record in SeqIO.parse(align_fa,"fasta"):
		if record.id==wt_name:
			continue
		rd[record.id]=str(record.seq[real_start:real_end])
	df_summary=pandas.DataFrame.from_dict(rd,orient='index',columns=["AA_seq"])
	#print(df_summary.shape)
	df_summary.loc["wt"]=[wt_seq[real_start:real_end]]
	print(df_summary)
	return df_summary
if __name__ == '__main__':
	#fl=glob.glob("results20200707/function/*_alleles_deal.xls")
	flagd={}
	#grnad={}
	#wtcdsd={}
	indeld={}
	locusd={}
	scored={}
	grnad={}
	for x in open("200825-FunctionalAnalysisInput-v6_deal.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		wt_cds=len(l[2])
		mut_cds=len(l[3])
		if wt_cds>mut_cds:
			flagd[l[0]]="delete"
		else:
			flagd[l[0]]="insert"
		#grnad[l[0]]=l[5].upper()
		#wtcdsd[l[0]]=l[2].upper()
		if l[-3]=="-":
			indeld[l[0]]=l[-3]
		else:
			indeld[l[0]]=float(l[-3])
		scored[l[0]]=float(l[-2])
		grnad[l[0]]=l[-1]
		locusd[l[0]]=l[1]
	extend=20
	fl=glob.glob("*_alleles_deal.xls")
	d=collections.defaultdict(dict)
	for fi in fl:
		print(fi)
		name=fi.split("/")[-1].split("_alleles_deal.xls")[0]
		print(name)
		sample=name.split("_")[1]
		#if name!="Rep1_54-051":
		#	continue
		df=pandas.read_csv(fi,sep="\t",header=0,index_col=0)
		d[name]["indel"]=indeld[sample]
		d[name]["Specificity_Score"]=scored[sample]
		d[name]["grna"]=grnad[sample]
		d[name]["locus"]=locusd[sample]
		p_sum=df[df["p_function"]!="-"]["#Reads"].sum()
		df_p=df[df['p_function']!="-"].copy()
		yes_sum=df[df["p_function"]=="yes"]["#Reads"].sum()
		all_sum=df["#Reads"].sum()
		sub_sum=df[(df["n_deleted"]==0)&(df["n_inserted"]==0)]["#Reads"].sum()
		indel_sum=all_sum-sub_sum
		mut_flag=flagd[sample]
		if mut_flag=="insert":
			inframe_sum=df[(df["n_inserted"]-df["n_deleted"]+1)%3==0]["#Reads"].sum()
		else:
			inframe_sum=df[(df["n_inserted"]-df["n_deleted"]-1)%3==0]["#Reads"].sum()
		if p_sum==0:
			d[name]["Functional_in_inframe"]=0
			d[name]["aa0"]=0
			d[name]["aa1"]=0
			d[name]["aa2"]=0
			d[name]["aa3"]=0
			d[name]["aa_others"]=0
		else:
			df_p['aa_change']=df_p["p_deleted"].astype(int)+df_p['p_inserted'].astype(int)+df_p['p_mutated'].astype(int)
			aa0=df_p[df_p["aa_change"]==0]["#Reads"].sum()
			aa1=df_p[df_p["aa_change"]==1]["#Reads"].sum()
			aa2=df_p[df_p["aa_change"]==2]["#Reads"].sum()
			aa3=df_p[df_p["aa_change"]==3]["#Reads"].sum()
			aa_others=df_p[df_p["aa_change"]>3]["#Reads"].sum()
			d[name]["Functional_in_inframe"]=yes_sum/float(indel_sum)
			d[name]["aa0"]=aa0/float(indel_sum)
			d[name]["aa1"]=aa1/float(indel_sum)
			d[name]["aa2"]=aa2/float(indel_sum)
			d[name]["aa3"]=aa3/float(indel_sum)
			d[name]["aa_others"]=aa_others/float(indel_sum)
		print(p_sum,inframe_sum,indel_sum)
		d[name]["inframe"]=p_sum/float(indel_sum)
		d[name]["inframe_old"]=inframe_sum/float(indel_sum)
		#print(d)
	df_r=pandas.DataFrame.from_dict(d,orient="index")
	df_r.to_csv("stat3.xls",sep="\t",index_label="sample")
	df_polar=df_r[df_r["indel"]!="-"].copy()
	df_polar.reset_index(drop=True,inplace=True)
	print(df_polar)
	df_polar.set_index(["locus","grna"],inplace=True)
	print(df_polar)
	df_polar['indel']=pandas.to_numeric(df_polar['indel'])
	df_polar_mean=df_polar.groupby(["grna","locus"]).mean()
	df_polar_mean["aaxd1"]=df_polar_mean["aa0"]+df_polar_mean["aa1"]
	df_polar_mean["aaxd2"]=df_polar_mean["aa0"]+df_polar_mean["aa1"]+df_polar_mean["aa2"]
	df_polar_mean["aaxd3"]=df_polar_mean["aa0"]+df_polar_mean["aa1"]+df_polar_mean["aa2"]+df_polar_mean["aa3"]
	print(df_polar_mean)
	print(df_polar_mean.columns)
	df_polar_mean.reset_index(drop=False,inplace=True)
	df_polar_mean["sort_grna"]=df_polar_mean["grna"].str.split("-",expand=True).iloc[:,3]
	df_polar_mean.sort_values(by=['sort_grna'],inplace=True)
	df_polar_r=df_polar_mean.loc[:,['locus','grna','inframe','Functional_in_inframe','aaxd1', 'aaxd2','aaxd3','indel','Specificity_Score','sort_grna']]
	df_polar_r.rename(columns={"inframe":'%Inframe',"Functional_in_inframe":"%Function_in_inframe","indel":"%Indel"},inplace=True)
	df_polar_r.loc[:,['%Inframe', '%Function_in_inframe', 'aaxd1', 'aaxd2','aaxd3']]=df_polar_r.loc[:,['%Inframe', '%Function_in_inframe', 'aaxd1', 'aaxd2','aaxd3']]*100
	df_polar_r.to_csv("polar3.xls",sep="\t",index=False)
