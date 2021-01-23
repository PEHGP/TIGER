#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import sys,re,os,pandas
import collections,re
import seaborn as sns
import matplotlib.pyplot as plt
import time,argparse
def get_his(inputfile,prefix,min_score):
	Fr=open("%s_HM_histogram.txt"%prefix,"w")
	Fr.write("indel_size\tCount\tHM_counts\n")
	dall=collections.defaultdict(int)
	dhm=collections.defaultdict(int)
	dr=pandas.DataFrame(columns=["indel_size","Count","HM_counts","HM_size"])
	for x in open(inputfile).readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		if l[4]==l[5]=="0":
			continue
		if l[9]=="True" and int(l[11])>min_score:
			dall[int(float(l[5]))-int(float(l[4]))]+=int(l[7])
			dhm[int(float(l[5]))-int(float(l[4]))]+=int(l[7])
			dd={'indel_size':int(float(l[5]))-int(float(l[4])),'Count':int(l[7]),'HM_counts':int(l[7]),'HM_size':len(l[10])}
		else:
			dall[int(float(l[5]))-int(float(l[4]))]+=int(l[7])
			dd={'indel_size':int(float(l[5]))-int(float(l[4])),'Count':int(l[7]),'HM_counts':0,'HM_size':0}
		dr=dr.append(dd,ignore_index=True)
	for x in sorted(dall.keys()):
		Fr.write(str(x)+"\t"+str(dall[x])+"\t"+str(dhm[x])+"\n")
	Fr.close()
	Fr2=open("%s_group_HM_size.txt"%prefix,"w")
	Fr2.write("indel_size\tHM_size\tCount\tHM_counts\n")
	for name,group in dr.groupby(['indel_size','HM_size']):
		Fr2.write("\t".join([str(n) for n in name])+"\t"+str(group['Count'].sum())+"\t"+str(group['HM_counts'].sum())+"\n")
	Fr2.close()
	dr.to_csv("%s_group_all.txt"%prefix,sep="\t",index=False)
def get_len_his(inputfile,outputfile,min_score):
	Fr=open(outputfile,"w")
	Fr.write("HM_size\tCount\n")
	d=collections.defaultdict(int)
	for x in open(inputfile).readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		if l[9]=='True' and int(l[11])>min_score:
			if l[4]==l[5]=="0":
				continue
			d[len(l[10])]+=int(l[7])
	for x in sorted(d.keys()):
		Fr.write(str(x)+"\t"+str(d[x])+"\n")
	Fr.close()
def plot_indel(indel_file,prefix):
	#bs=[-30,-20]+range(-10,11)+[20,30]
	bs=range(-40,11)
	df=pandas.read_csv(indel_file,sep="\t",header=0,index_col=None)
	if df["Count"].sum():
		df['Count_freq']=df["Count"]/df["Count"].sum()
		df['HM_counts_freq']=df["HM_counts"]/df["Count"].sum()
	else:
		df.loc[0]=[0,0,0]
		df.loc[:,"Count_freq"]=0
		df.loc[:,"HM_counts_freq"]=0
	df2=pandas.DataFrame()
	for i,bi in enumerate(bs):
		if i==0:
			label="<"+str(bi)
			df2[label]=df[df['indel_size']<bi].sum()
			label="[%s,%s)"%(bs[i],bs[i+1])
			df2[label]=df[(df['indel_size']>=bs[i])&(df['indel_size']<bs[i+1])].sum()
		elif i==len(bs)-1:
			label=">="+str(bi)
			df2[label]=df[df['indel_size']>=bi].sum()
		else:
			label="[%s,%s)"%(bs[i],bs[i+1])
			df2[label]=df[(df['indel_size']>=bs[i])&(df['indel_size']<bs[i+1])].sum()
	df2=df2.T.reset_index()
	f, ax = plt.subplots(2,1,figsize=(17,10))
	sns.set_context("paper")
	sns.set_style("white")
	sns.set_color_codes("pastel")
	a1=sns.barplot(x="indel_size", y="Count_freq",data=df, color="b",label='Homology-less',ax=ax[0])
	sns.barplot(x="indel_size", y="HM_counts_freq",data=df, color="r",label='MMEJ',ax=ax[0])
	ax[0].legend(ncol=1, loc="upper right", frameon=True)
	if df["Count"].sum()!=0:
		ax[0].text(0,df["Count_freq"].max()-0.1, u'Indelsum=%s\nHMsum=%s\nHM(%%)=%s%%'%(df["Count"].sum(),df["HM_counts"].sum(),float(df["HM_counts"].sum())/df["Count"].sum()*100))
	ax[0].set(ylabel="Frequency(%)",xlabel="Indel size")
	#print dir(a1)
	"""
	for p in a1.patches:
		print p.get_height()
		ax[0].text(p.get_x() + p.get_width()/2.,p.get_height(),p.get_x(), fontsize=7,color='red',ha='center',va='bottom')
	"""
	#-------------------------------------------------------------------
	sns.barplot(x="index", y="Count_freq",data=df2, color="b",label='Homology-less',ax=ax[1])
	sns.barplot(x="index", y="HM_counts_freq",data=df2, color="r",label='MMEJ',ax=ax[1])
	ax[1].legend(ncol=1, loc="upper right", frameon=True)
	if df2["Count"].sum()!=0:
		ax[1].text(0,df2["Count_freq"].max()-0.1, u'Indelsum=%s\nHMsum=%s\nHM(%%)=%s%%'%(df2["Count"].sum(),df2["HM_counts"].sum(),float(df2["HM_counts"].sum())/df2["Count"].sum()*100))
	ax[1].set(ylabel="Frequency(%)",xlabel="Indel size")
	f.savefig("%s_indel_size.png"%prefix,format="png",dpi=300)
	plt.close(f)
	df2.to_csv("%s_plot_indel.txt"%prefix,sep="\t",index=False)
def plot_size(hmsize_file,prefix):
	df=pandas.read_csv(hmsize_file,sep="\t",header=0,index_col=None)
	if df['Count'].sum():
		df['Count_freq']=df['Count']/df['Count'].sum()
	else:
		df.loc[0]=[0,0]
		df.loc[:,"Count_freq"]=0
	f, ax = plt.subplots()
	sns.set_context("paper")
	sns.set_style("white")
	sns.barplot(x="HM_size", y="Count_freq",data=df, color="b")
	ax.set(ylabel="HM_dependen_deletion(%)",xlabel="Indel size")
	f.savefig("%s_hm_size.png"%prefix,format="png",dpi=300)
	plt.close(f)
def findhomology(ref,start,end):
	ref=ref.upper()
	t=0
	homoseq1=""
	while True:
		#print ref
		#print start
		#print end
		if start+t<end and ref[start+t]==ref[end+t]:
			homoseq1+=ref[start+t]
			t+=1
		else:
			break
	t=1
	homoseq2=""
	while True:
		if end-t>=start and ref[start-t]==ref[end-t]:
			homoseq2+=ref[start-t]
			t+=1
		else:
			break
	return homoseq1,homoseq2[::-1]
def judge(flag,line_list):
	#print flag
	if flag=="default":
		return line_list[3]=="MODIFIED"
	else:
		return line_list[3]=="MODIFIED" and (line_list[2]==flag or line_list[2]=="AMBIGUOUS")
	#else:
	#	print "=============flag not right:%s=============="%flag
	#	print "error error error error error"
	#	sys.exit()
def dealseq(seq):
	dseq=re.sub('[NATCG]','-',seq,10)
	#print dseq
	dseq=re.sub('[NATCG]','-',dseq[::-1],10)
	return dseq[::-1]
def get_region(line,grna,ptime,extend): #list,Reverse complementarity yes
	if not grna:
		print("====grna is not exist!====")
		sys.exit()
	Fr=open("temp"+ptime+".fa","w")
	Fr.write(">temp"+ptime+"\n")
	Fr.write(line[1]+"\n")
	Fr.close()
	Fr=open("grna"+ptime+".fa","w")
	Fr.write(">grna"+ptime+"\n")
	Fr.write(grna+"\n")
	Fr.close()
	os.system("bwa index temp%s.fa"%ptime)
	start=-1
	k=12
	while start<0:
		if k<5:
			break
		lines=os.popen("bwa mem -t 30 -k %s -T 0 temp%s.fa grna%s.fa"%(k,ptime,ptime)).readlines()
		l=lines[2].rstrip().split("\t")
		start=int(l[3])-1
		end=start+len(grna)
		k=k-1
	#print line[1][start:end]
	#print grna
	#print line[0][start-extend:end+extend]
	if start<0:
		print("====not found grna region====")
		print(line)
		print(grna)
		sys.exit()
	return start-extend,end+extend
def main(f,prefix,ref,min_score):#deletion judge
	ptime=str(time.time()).split(".")[0]
	Fr1=open("%s_Alleles_frequency_table_HM.txt"%prefix,"w")
	Fr2=open("%s_homology_special.txt"%prefix,"w")
	Fr3=open("%s_log_HM.txt"%prefix,"w")
	lines=open(f).readlines()
	Fr1.write(lines[0].rstrip()+"\thomology\thomoseq\thomoseq_score\tall_info\tline_num\n")
	scored={"A":2,"T":2,"G":3,"C":3,"N":0}
	extend=20
	grna=re.search("-g (.*?) ",open("%s_CRISPResso_RUNNING_LOG.txt"%prefix).read()).group(1)
	print(grna)
	for n,x in enumerate(lines):
		x=x.rstrip()
		l=x.split("\t")
		#if l[3]=="MODIFIED":
		if judge(ref,l):
			if len(l[0])==len(l[1]):
				r=[]
				if int(l[4])>0:
					gr=get_region(l,grna,ptime,extend)
					for m in re.finditer("-+",l[0]): #need change
						Fr3.write(str(n+1)+"\t"+str(m.start())+"\t"+str(m.end())+"\t"+str(m.group())+"\t"+str(len(m.group()))+"\n") #m.end()-1 is the last "-"
						#print l[0]
						#print m.group(),m.start(),m.end(),n
						if m.start()==0 or m.end()==len(l[0]):
							continue
						#print l[0][m.start()-1],l[0][m.end()]
						#l[1][m.start()]
						#l[1][m.end()]
						jo=max(0,min(m.end(), gr[1]) - max(m.start(),gr[0]))
						if jo==0:
							continue
						h1,h2=findhomology(l[1],m.start(),m.end())
						s1=sum([scored[s] for s in h1])
						s2=sum([scored[s] for s in h2])
						#print h1,h2
						if s1!=0 or s2!=0:
								r.append((h1,h2,s1,s2,len(m.group())))
				if len(r)==1:
					if r[0][2]>=r[0][3]:
						if len(r[0][0])<=r[0][4]:
							flag=True
						else:
							flag=False
						Fr1.write(x+"\t%s\t"%flag+r[0][0]+"\t"+str(r[0][2])+"\t"+",".join([str(ri) for ri in r])+"\t"+str(n)+"\n")
					elif r[0][2]<r[0][3]:
						if len(r[0][1])<=r[0][4]:
							flag=True
						else:
							flag=False
						Fr1.write(x+"\t%s\t"%flag+r[0][1]+"\t"+str(r[0][3])+"\t"+",".join([str(ri) for ri in r])+"\t"+str(n)+"\n")
					#else:
					#	Fr1.write(x+"\tTrue\t-\t%s\t"%(r[0][1])+",".join([str(ri) for ri in r])+"\n") #False?
				elif len(r)>1:
					#fm=""
					r2=""
					score=0
					gaplen=0
					for ri in r:
						if ri[2]>score:
							r2=ri[0]
							score=ri[2]
							gaplen=ri[4]
						if ri[3]>score:
							r2=ri[1]
							score=ri[3]
							gaplen=ri[4]
						#fm+=",".join([str(rii) for rii in ri])+"||"

					if score>0 and len(r2)<=gaplen:
						Fr1.write(x+"\tTrue\t%s\t%s\t"%(r2,score)+",".join([str(ri) for ri in r])+"\t"+str(n)+"\n")
					else:
						Fr1.write(x+"\tFalse\t-\t0\t%s\t%s\n"%(",".join([str(ri) for ri in r]),n))
				else:
					Fr1.write(x+"\tFalse\t-\t0\t-\t%s\n"%(n))
			else:
				Fr2.write(x+"length not same:%s\n"%(n))
				Fr3.write("%s,length not same:%s,%s\n"%(n,len(l[0]),len(l[1])))
	Fr3.close()
	Fr2.close()
	Fr1.close()
	get_his("%s_Alleles_frequency_table_HM.txt"%prefix,prefix,min_score)
	get_len_his("%s_Alleles_frequency_table_HM.txt"%prefix,"%s_HM_size_histogram.txt"%prefix,min_score)
	plot_indel("%s_HM_histogram.txt"%prefix,prefix)
	plot_size("%s_HM_size_histogram.txt"%prefix,prefix)
def get_parser():
	parser = argparse.ArgumentParser(description="deal2_new20190723.py",formatter_class=argparse.RawTextHelpFormatter,epilog='hehe')
	parser.add_argument('--prefix',dest="prefix",type=str,help="results prefix")
	parser.add_argument('--ref',dest="ref", type=str,default='default',help="ref name")
	parser.add_argument('--min-score',dest="score",type=int,default=3,help="min homology score  default:3")
	parser.add_argument('--debug',dest="debug",action="store_true",help="debug program")
	return parser
if __name__=="__main__":
	if len(sys.argv)==1:
		get_parser().print_help()
		sys.exit()
	else:
		args = get_parser().parse_args()
	if os.path.isdir(args.prefix):
		print("%s is exist"%args.prefix)
		sys.exit()
	os.system("mkdir %s"%args.prefix)
	samplelist=[]
	lines=os.popen("find -L . -name \"Alleles_frequency_table.zip\"").readlines()
	for xi in lines:
		pa=xi.rstrip()
		pl=pa.split("/")
		#print pl 
		#sys.exit()
		print(pl)
		if len(pl)==4: #need change
			print(pl[2])
			sample_name=pl[2].split("_")[-1]
			if not sample_name in samplelist:
				samplelist.append(sample_name)
			else:
				print(sample_name)
				print("========same sample name=======")
				print("error error error error error")
				sys.exit()
			print(pa)
			#sys.exit()
			os.system("unzip -o %s -d %s"%(pa,args.prefix))
			os.system("mv %s/Alleles_frequency_table.txt %s/%s_Alleles_frequency_table.txt"%(args.prefix,args.prefix,sample_name))
			os.system("cp %s/CRISPResso_RUNNING_LOG.txt %s/%s_CRISPResso_RUNNING_LOG.txt"%("/".join(pl[:-1]),args.prefix,sample_name))
	path=os.getcwd()
	os.chdir(path+"/"+args.prefix)
	for s in samplelist:
		print(s)
		main("%s_Alleles_frequency_table.txt"%s,s,args.ref,int(args.score))
		if not args.debug:
			os.remove("%s_Alleles_frequency_table.txt"%s)
			os.remove("%s_homology_special.txt"%s)
			os.remove("%s_log_HM.txt"%s)
			os.remove("%s_CRISPResso_RUNNING_LOG.txt"%s)
			os.system("rm -rf temp*")
			os.system("rm -rf grna*.fa")
	os.chdir(path)
