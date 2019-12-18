#!/usr/bin/env python2.7
import pandas,sys
f=sys.argv[1] #indel_size_freq.txt
df=pandas.read_table(f,sep="\t",header=0,index_col=0)
df2=df[((df.index+1)%3==0)].copy() #+1 need change
df2.loc["inframe",:]=df2.sum()
il=[5,2,-1,-4,-7] #need change
df2.loc["others",:]=df2[~df2.index.isin(il+["inframe"])].sum() 
df2.loc[["inframe"]+il+["others"],:].to_csv("inframe.txt",sep="\t")
