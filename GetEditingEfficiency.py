#!/usr/bin/env python
import sys,os,pandas
line=os.popen("find . -name 'Alleles_frequency_table.zip' ").readlines()
d={}
for x in line:
	x=x.rstrip()
	#os.system("unzip -o %s"%x)
	df=pandas.read_csv(x,sep="\t",header=0,index_col=None)
	alreads=df["#Reads"].sum()
	modified=df[df["Read_Status"]=="MODIFIED"]["#Reads"].sum()
	sub=df[(df["Read_Status"]=="MODIFIED")&(df["n_deleted"]==0)&(df["n_inserted"]==0)]["#Reads"].sum()
	d[x]=[alreads,modified,sub]
df2=pandas.DataFrame.from_dict(data=d,orient='index',columns=["Reads_aligned","Modified","Only_substitutions"])
df2["indel_counts"]=df2["Modified"]-df2["Only_substitutions"]
df2["indel_freq"]=(df2["Modified"]-df2["Only_substitutions"])/df2["Reads_aligned"]*100
df2.to_csv("EditingEfficiency.txt",sep="\t",index=True,index_label="sample")
