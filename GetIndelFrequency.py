#!/usr/bin/env python
import os,pandas
def get_indel(table):
	df=pandas.read_csv(table,sep="\t",header=0,index_col=None)
	df2=df[(df["n_deleted"]!=0)|(df["n_inserted"]!=0)]
	df2['indel']=df2['n_inserted']-df2['n_deleted']
	r=df2.groupby('indel')["#Reads"].sum()
	r=r/r.sum()*100
	return r
def batch_indel(sample,grna_name):
	j3=["53-051","54-051","55-051","56-051"]
	sl=os.popen("find -L . -name \"CRISPResso_on_%s\""%(sample)).readlines()
	rl=[]
	for p in sl:
		p=p.rstrip()
		if sample in sl:
			pa=p+"/Pure_Alleles_frequency_table.txt"
			r=get_indel(pa)
		else:
			pa=p+"/Alleles_frequency_table.zip"
			os.system("unzip -o %s"%(pa))
			r=get_indel("Alleles_frequency_table.txt")
		rl.append(r)
	df=pandas.concat(rl,axis=1,join='outer').fillna(0)
	#print df
	df=df.mean(axis=1)
	df=df.to_frame()
	#print df
	df.columns=["Predicted frequency"]
	#print df
	df.to_csv("%s_%s_indel.xls"%(sample,grna_name),sep="\t",index_label="Indel length")
if __name__ == '__main__':
	for x in open("200709-inDelphi-input.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		#l[4] 001 l[5] 051
		batch_indel(l[4],l[1])
		batch_indel(l[5],l[1])
		#break
