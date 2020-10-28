#!/usr/bin/env python
import sys,os
import pandas,re
df_list=[]
freq_list=[]
sample_list=[]
df2_list=[]
freq2_list=[]
labels=["Count","HM_counts","HM_less"]
lines=os.popen("find . -name \"CRISPResso_on_*_HM_histogram.txt\"").readlines()
#print lines
if not lines:
	print x
	print "error error error error"
	sys.exit()
for i,f in enumerate(lines):
	f=f.rstrip()
	#print f
	pr=f.split("/")[1].replace("_Output","")
	m=re.search("CRISPResso_on_(.*?)_HM_histogram",f)
	name=m.group(1)
	f2="/".join(f.split("/")[:-1])+"/"+"CRISPResso_on_%s_HM_size_histogram.txt"%name
	print name,pr
	print f
	print f2
	xi=pr+"_"+name
	print xi
	#sys.exit()
	sample_list.append(xi)
	df=pandas.read_csv(f,sep="\t",header=None,index_col=0,skiprows=1,names=["%s_Count"%xi,"%s_HM_counts"%xi])
	df2=pandas.read_csv(f2,sep="\t",header=None,index_col=0,skiprows=1,names=["%s_Count"%xi])
	#print df[xi+"_Count"]
	df_freq=pandas.DataFrame()
	df2_freq=pandas.DataFrame()
	df[xi+"_HM_less"]=df[xi+"_Count"]-df[xi+"_HM_counts"]
	s=df[xi+"_Count"].sum()
	df2_freq[xi+"_Count"]=df2[xi+"_Count"]/s
	for la in labels:
		df_freq[xi+"_"+la]=df[xi+"_"+la]/s
	df_list.append(df)
	freq_list.append(df_freq)
	df2_list.append(df2)
	freq2_list.append(df2_freq)
	#break
print len(set(sample_list))
df_r=pandas.concat(df_list,axis=1)
df_r=df_r.fillna(0)
freq_r=pandas.concat(freq_list,axis=1)
freq_r=freq_r.fillna(0)
df2_r=pandas.concat(df2_list,axis=1)
df2_r=df2_r.fillna(0)
freq2_r=pandas.concat(freq2_list,axis=1)
freq2_r=freq2_r.fillna(0)
df2_r.to_csv("HM_size_histogram.txt",sep="\t",index_label='indel_size')
freq2_r.to_csv("HM_size_histogram_freq.txt",sep="\t",index_label='indel_size')
for la in labels:
	sl=[]
	for s in sample_list:
		sl.append(s+"_"+la)
	df_r[sl].to_csv("HM_histogram_%s.txt"%la,sep="\t",index_label="indel_size")
	freq_r[sl].to_csv("HM_histogram_freq_%s.txt"%la,sep="\t",index_label="indel_size")
