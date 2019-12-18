#!/usr/bin/env python2.7
import sys,os,pandas
def GetIndelOne(f,name):
	df=pandas.read_csv(f,sep="\t",header=0,index_col=None)
	df["indel_size"]=df["n_inserted"]-df["n_deleted"]
	sub=df[(df["Read_Status"]=="MODIFIED")&(df["n_deleted"]==0)&(df["n_inserted"]==0)]["#Reads"].sum()
	modified=df[df["Read_Status"]=="MODIFIED"]["#Reads"].sum()
	indel_count=modified-sub
	df2=df[(df["Read_Status"]=="MODIFIED")&((df["n_deleted"]!=0)|(df["n_inserted"]!=0))]
	df3=df2.groupby('indel_size')['#Reads'].sum().to_frame()
	df3.columns=[name]
	df3_freq=(df2.groupby('indel_size')['#Reads'].sum()/indel_count*100).to_frame()
	df3_freq.columns=[name+"_freq"]
	return df3,df3_freq
def GetIndelList(flist):
	rl=[]
	rl_freq=[]
	for x in flist:
		x=x.rstrip()
		name=os.path.dirname(x).replace("./","").replace("/","_")
		df3,df3_freq=GetIndelOne(x,name)
		rl.append(df3)
		rl_freq.append(df3_freq)
	df_r=pandas.concat(rl,axis=1)
	df_r=df_r.fillna(0)
	df_freq=pandas.concat(rl_freq,axis=1)
	df_freq=df_freq.fillna(0)
	return df_r,df_freq

if __name__ == '__main__':
	line=os.popen("find . -name 'Pure_Alleles_frequency_table.txt' ").readlines()
	df_r,df_freq=GetIndelList(line)
	df_r.to_csv("indel_size.txt",sep="\t")
	df_freq.to_csv("indel_size_freq.txt",sep="\t")

