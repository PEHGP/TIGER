#!/usr/bin/env python
import sys,os,pandas
from Bio import SeqIO
import re,t7
#from CRISPResso2 import CRISPRessoShared
from Bio.Seq import Seq
import glob,collections
from multiprocessing import Pool as ProcessPool
def final_pip(params):
	label,name,wt_cds_seq,mutant_cds_seq,grna_seq,target_aa,inframed,df2=params
	rep=os.popen("awk '{print $8}' %s_blast.txt|sort|uniq -c|sort -n -r"%name).readlines()[0].rstrip().split()[1]
	rl=os.popen("awk '{print $6}' %s_blast.txt|sort -u"%name).readlines()[0].rstrip()
	if int(rep)!=int(rl):
		print("algn error please check %s_blast.txt"%name)
		sys.exit()
	rsp=os.popen("awk -F'\t' '{print $7}' %s_blast.txt|sort|uniq -c|sort -k 1 -n -r"%name).readlines()[0].rstrip().split()[1]
	print("ref start align position: %s"%rsp)
	t7.deal_blast(name,rep,rsp,wt_cds_seq)
	#os.system("mafft --globalpair --thread 60 --inputorder %s_deal.fa>%s_temp_mafft.fa"%(name,name))
	grna_seq=t7.get_grna(grna_seq,label)
	cut_point,orientation=t7.my_cut(grna_seq,wt_cds_seq)
	Fr_aa=open(name+".faa","w")
	Fr_log=open(name+"_inframeerror.fa","w")
	Fr_aa.write(">wt_cds\n")
	Fr_aa.write(str(Seq(wt_cds_seq).translate())+"\n")
	for record in SeqIO.parse("%s_deal.fa"%name,"fasta"):
		Fr=open("%s_single_temp.fa"%name,"w")
		Fr.write(">wt_cds\n")
		Fr.write(wt_cds_seq+"\n")
		Fr.write(">mutant_cds\n")
		Fr.write(mutant_cds_seq+"\n")
		Fr.write(">"+record.id+"\n")
		Fr.write(str(record.seq)+"\n")
		Fr.close()
		os.system("mafft --quiet %s_single_temp.fa>%s_single_temp_mafft.fa"%(name,name))
		if not t7.judge_mut(name):
			os.system("mafft --quiet --genafpair --thread 3 --inputorder %s_single_temp.fa>%s_single_temp_mafft.fa"%(name,name))
		if not t7.judge_mut(name):
			os.system("mafft --op 0 --ep 1 --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
		if not t7.judge_mut(name):
			os.system("mafft --op 0 --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
		if not t7.judge_mut(name):
			os.system("mafft --globalpair --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
		mafft_handle=SeqIO.index("%s_single_temp_mafft.fa"%(name),"fasta")
		reald,mut_flag=t7.mut_reads(cut_point,mafft_handle)
		t7.translate_and_align(reald,inframed,mut_flag,name,Fr_aa,Fr_log)
	Fr_aa.close()
	Fr_log.close()
	faa_num=int(os.popen("wc -l %s.faa"%name).readlines()[0].rstrip().split()[0])
	if faa_num==2:
		os.system("cp %s.faa %s_align.faa"%(name,name))
	else:
		os.system("mafft --quiet --thread 40 --globalpair %s.faa >%s_align.faa"%(name,name))
	aa_file="mut_cds.fa"
	#target_aa=l[-1].upper()
	extend=20
	grna_prot=t7.get_grna_prot(wt_cds_seq,grna_seq,extend,cut_point)
	consurf_file=t7.get_consurf_file(target_aa,aa_file)
	#print target_aa
	t7.get_var(df2,"%s_align.faa"%name,grna_prot,consurf_file,name)
if __name__ == '__main__':
	finish_l=[]
	paramsl=[]
	for ri in glob.glob("*_alleles_deal.xls"):
		finish_l.append(ri.split("_alleles_deal.xls")[0])
	for x in open("Input-051.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		grna_seq=l[5].upper()
		label=l[0]
		wt_cds_seq=l[2].upper()
		mutant_cds_seq=l[3].upper()
		target_aa=l[-1].upper()
		sl=os.popen("find -L . -name \"CRISPResso_on_%s\""%(l[0])).readlines()
		print(sl)
		for p in sl:
			p=p.rstrip()
			name=p.split("/")[1].split("-")[0]+"_"+l[0]
			if name in finish_l:
				continue
			Fr_cds=open(name+"_wtcds.fa","w")
			Fr_cds.write(">wt_cds\n")
			Fr_cds.write(l[2].upper()+"\n")
			Fr_cds.close()
			os.system("makeblastdb -in %s_wtcds.fa -out %s_wtcds -dbtype nucl"%(name,name))
			#Fr.write(">mutant_cds\n")
			#Fr.write(l[3].upper()+"\n")
			pa=p+"/Alleles_frequency_table.zip"
			print(pa)
			#os.system("unzip -o %s"%(pa))
			df=pandas.read_csv(pa,sep="\t",header=0,index_col=None,compression='zip')
			df2=df[(df["Read_Status"]=="MODIFIED")]
			Fr=open(name+"_temp.fa","w")
			inframed={}
			for index,row in df2.iterrows():
				seq=row["Aligned_Sequence"].replace("-","")
				seq_name="temp_%s"%(index)
				Fr.write(">"+seq_name+"\n")
				Fr.write(seq+"\n")
				inframed[seq_name]=row["n_inserted"]-row["n_deleted"]
			Fr.close()
			os.system("blastn -db %s_wtcds -query %s_temp.fa -out %s_blast.txt -num_threads 30 -evalue 1e-10 -outfmt \"6 qseqid qlen qstart qend sseqid slen sstart send qcovs bitscore evalue pident\""%(name,name,name))
			paramsl.append((label,name,wt_cds_seq,mutant_cds_seq,grna_seq,target_aa,inframed,df2))
			
			#break
		#break
	pool = ProcessPool(10)
	pool.map(final_pip,paramsl)
	pool.close()
	pool.join()
