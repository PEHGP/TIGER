#!/usr/bin/env python
import sys,os,pandas
import collections
if __name__ == '__main__':
	line=[]
	line=os.popen("find . -name 'Alleles_frequency_table_around_sgRNA*.txt'").readlines()
	#line=os.popen("find . -name 'Pure_Alleles_frequency_table_around_sgRNA*.txt' ").readlines()
	for x in line:
		x=x.rstrip()
		seql=[]
		dfl=[]
		print x
		df=pandas.read_csv(x,sep="\t",header=0,index_col=None)
		df["indel_size"]=df["n_inserted"]-df["n_deleted"]
		#df["inframe"]=(df["indel_size"]+1)%3==0
		df["inframe"]="hehe"
		modified=df[df["Unedited"]==False]["#Reads"].sum()
		sub=df[(df["Unedited"]==False)&(df["n_deleted"]==0)&(df["n_inserted"]==0)]["#Reads"].sum()
		indel_count=modified-sub
		df["freq_total"]=df["#Reads"]/indel_count*100
		df2=df[(df["Unedited"]==False)&((df["n_deleted"]!=0)|(df["n_inserted"]!=0))].copy()
		for name,dfg in df2.groupby("indel_size"):
			dfg=dfg.copy()
			dfg["freq_group"]=dfg["#Reads"]/dfg["#Reads"].sum()*100
			#print name,dfg["#Reads"].sum()
			#print dfg["reads_group"]
			#print score_total
			score_total=""
			score_group=""
			d=collections.defaultdict(int)
			dfl.append(dfg)
			for index,row in dfg.iterrows():
				#print row
				align_seq=list(row["Aligned_Sequence"])
				ref_seq=list(row["Reference_Sequence"])
				seql.append(str(row["indel_size"])+"\tseq\t"+"\t".join(align_seq)+"\t"+str(row["freq_total"]))
				for i,s in enumerate(align_seq):
					if name<0 and s=="-":
						d[i]+=row["#Reads"]
					if name>0 and ref_seq[i]=="-":
						d[i]+=row["#Reads"]
			#print score_total
			#print score_group
			#print d
			score_total+=str(row["indel_size"])+"\tscore_total\t"
			score_group+=str(row["indel_size"])+"\tscore_group\t"
			for i in range(len(row["Aligned_Sequence"])):
				score_group+=str(float(d[i])/dfg["#Reads"].sum()*100)+"\t"
				score_total+=str(float(d[i])/indel_count*100)+"\t"
			#print score_group
			seql.append(score_total+"None")
			seql.append(score_group+"None")
		name=x.replace("./","").replace("/","_").split("sgRNA")[0].replace("_Alleles_frequency_table_around_","")
		if dfl:
			df3=pandas.concat(dfl)
			df3.loc[:,["Aligned_Sequence","Reference_Sequence","Unedited","n_deleted","n_inserted","n_mutated","#Reads","freq_total","freq_group","indel_size","inframe"]].to_csv(name+"_Editing_product.txt",sep="\t",index=False)
			Fr=open(name+"_Editing_product_seqlist.txt",'w')
			Fr.write("\n".join(seql)+"\n")
			Fr.close()
		else:
			os.system("touch %s_Editing_product.txt"%name)
			os.system("touch %s_Editing_product_seqlist.txt"%name)
		#break
