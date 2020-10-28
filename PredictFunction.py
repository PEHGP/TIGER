#!/usr/bin/env python
import sys,os,pandas
from Bio import SeqIO
import re
#from CRISPResso2 import CRISPRessoShared
from Bio.Seq import Seq
import glob,collections
def get_region(cut_point,mafft_handle):
	wt_seq=str(mafft_handle["wt_cds"].seq)
	point=[]
	for i,c in enumerate(wt_seq):
		if c!="-":
			point.append(i)
	mafft_cut_point=point[cut_point]
	start=point[0]
	#start=0
	end=point[-1]+1
	#end=len(wt_seq)
	return start,end,mafft_cut_point
def mut_reads(cut_point,mafft_handle):
	start,end,mafft_cut_point=get_region(cut_point,mafft_handle)
	wt_seq=str(mafft_handle["wt_cds"].seq)
	mut_seq=str(mafft_handle["mutant_cds"].seq)
	pm=[]
	t=0
	for w,m in zip(wt_seq,mut_seq):
		if w!=m:
			pm.append(t)
			if w=="-":
				mut_flag="insert"
			elif m=="-":
				mut_flag="delete"
			else:
				print("mut error")
				print(t,w,m)
				print(wt_seq)
				print(mut_seq)
				sys.exit()
		t+=1
	if len(pm)>1:
		print("========mut is not one=======")
		sys.exit()
	reald={}
	for rn in mafft_handle:
		if rn=="wt_cds":
			reald[rn]=list(mafft_handle[rn].seq)[start:end]
			continue
		if rn=="mutant_cds":
			reald[rn]=list(mafft_handle[rn].seq)[start:end]
			continue
		p_mut=pm[0]
		reads_seq=list(mafft_handle[rn].seq)
		if mut_flag=="delete":
			if reads_seq[p_mut]!="-": #need change?
				reads_seq[p_mut]="-"
			elif reads_seq[p_mut]=="-":
				mut_start=p_mut
				mut_end=p_mut
				while reads_seq[mut_start]=="-":
					mut_start=mut_start-1
					if mut_start<0:
						mut_start=0
						break
				while reads_seq[mut_end]=="-":
					mut_end=mut_end+1
					if mut_end==len(reads_seq):
						mut_end=mut_end-1
						break
				if mafft_cut_point>p_mut:
					if reads_seq[mut_end]!="-":
						reads_seq[mut_end]="-"
					else:
						reads_seq[start_end]="-"
				elif mafft_cut_point<=p_mut: #care =
					if reads_seq[mut_start]!="-":
						reads_seq[mut_start]="-"
					else:
						reads_seq[mut_end]="-"
			real_seq=reads_seq[start:end]
		elif mut_flag=="insert":
			if reads_seq[p_mut]=="-":
				mut_start=p_mut
				mut_end=p_mut
				while reads_seq[mut_start]=="-":
					mut_start=mut_start-1
					if mut_start<0:
						mut_start=0
						break
				while reads_seq[mut_end]=="-":
					mut_end=mut_end+1
					if mut_end==len(reads_seq):
						mut_end=mut_end-1
						break
				if mut_start==mut_end==p_mut:
					reads_seq[p_mut]=mut_seq[p_mut]
				else:
					if abs(mafft_cut_point-mut_start)>=abs(mafft_cut_point-mut_end):
						reads_seq[mut_start+1]=mut_seq[mut_start+1]
					else:
						reads_seq[mut_end-1]=mut_seq[mut_end-1]
				real_seq=reads_seq[start:end]
			elif reads_seq[p_mut]!="-":
				if mafft_cut_point>p_mut:
					real_seq=reads_seq[start:p_mut+1]+[mut_seq[p_mut]]+reads_seq[p_mut+1:end]
				elif mafft_cut_point<=p_mut: #care =
					real_seq=reads_seq[start:p_mut]+[mut_seq[p_mut]]+reads_seq[p_mut:end]
		reald[rn]=real_seq
	return reald,mut_flag
def translate_and_align(reald,inframed,mut_flag,preifx,Fr,Fr_log):
	#pl=[]
	for rn in reald:
		if rn=="wt_cds" or rn=="mutant_cds":
			continue
		seq=("".join(reald[rn])).replace("-","")
		if len(seq)%3==0:
			if mut_flag=="insert" and (inframed[rn]+1)%3!=0:
				#print "inframe error"
				Fr_log.write(">"+rn+"\n")
				Fr_log.write("".join(reald[rn])+"\n")
			if mut_flag=="delete" and (inframed[rn]-1)%3!=0:
				#print "inframe error"
				Fr_log.write(">"+rn+"\n")
				Fr_log.write("".join(reald[rn])+"\n")
			#pl.append(rn)
			Fr.write(">"+rn+"\n")
			Fr.write(str(Seq(seq).translate())+"\n")
def get_grna_prot(wt_cds,grna_seq,extend,cut_point):
	win=3
	cdsd={}
	for start,end in [(i,i+win) for i in range(0,len(wt_cds),win)]:
		#print start,end-1
		for ti in range(start,end):
			cdsd[ti]=(start,end)
	start,end=cut_point-extend,cut_point+extend
	if start<0:
		start=0
	if end>len(wt_cds):
		end=len(wt_cds)
	seq_deal=Seq(wt_cds[cdsd[start][0]:cdsd[end-1][1]])
	grna_prot=seq_deal.translate()
	return str(grna_prot).upper()
def get_consurf(consurf,wt_seq):
	conl=[]
	seq=""
	for x in open(consurf):
		x=x.rstrip()
		l=x.split("\t")
		conl.append(l)
		seq=seq+l[0]
	m=re.search(wt_seq,seq)
	start,end=m.span()
	return conl[start:end]
def get_consurf_file(target_aa,aa_file):
	aa_name=""
	for record in SeqIO.parse(aa_file,"fasta"):
		if str(record.seq.upper())==target_aa:
			aa_name=record.id
			break
	consurf_file="consurf/%s_consurf_deal.txt"%aa_name
	return consurf_file

def get_var(df_modified,align_fa,grna_prot,consurf_file,prefix):
	acidic_aa=["D","E"]
	basic_aa=["R","K","H"]
	wt_name="wt_cds"
	seqd=SeqIO.index(align_fa,"fasta")
	wt_seq=str(seqd[wt_name].seq).upper()
	wt_seq_deal=wt_seq.replace("-","").upper()
	conl=get_consurf(consurf_file,wt_seq_deal) #need change
	m=re.search(grna_prot,wt_seq_deal)
	if m:
		start,end=m.span() #end is not include
	else:
		print("grna_prot match error")
		sys.exit()
	pseq=[]
	conl_align=[]
	t=0
	for n,s in enumerate(wt_seq):
		if s!="-":
			pseq.append(n)
			conl_align.append(conl[t])
			t+=1
		else:
			conl_align.append(["-",0,"0"])
	#wt_seq_deal=wt_seq.replace("-","").upper()
	real_start=pseq[start]
	if end==len(pseq):
		real_end=pseq[end-1]+1
	else:
		real_end=pseq[end]
	rd={}
	for record in SeqIO.parse(align_fa,"fasta"):
		if record.id==wt_name:
			continue
		delete=0
		substitute=0
		insert=0
		t=0
		fuc="yes"
		con_num=0
		con_ch=0
		uncon_num=0
		uncon_ch=0
		print(conl_align[real_start:real_end])
		print(real_start,real_end,len(conl_align))
		for w,v,sr in zip(wt_seq[real_start:real_end],record.seq[real_start:real_end],conl_align[real_start:real_end]):
			print(t,w,v,sr)
			saa=sr[2].split(",")
			#clas=get_class(saa)
			if int(sr[1])>7:
				con_num+=1
			else:
				uncon_num+=1
			if w!=v:
				if int(sr[1])>7:
					con_ch+=1
				else:
					uncon_ch+=1
			if v=="*":
				fuc="no"
			if w!=v and w!="-" and v!="-":
				substitute+=1
				if int(sr[1])>7 and (not v in saa):
					if (w in acidic_aa) and (v in acidic_aa):
						pass
					elif (w in basic_aa) and (v in basic_aa):
						pass
					else:
						fuc="no"
				#else:
					#if (not v in saa) and (not v in clas):
				#	if not v in saa:
				#		fuc="no"
			elif w!=v and w=='-' and v!="-":
				insert+=1
				it=t
				while it>=0 and conl_align[real_start:real_end][it][0]=="-":
					it=it-1
				if it<0:
					b=8
				else:
					b=conl_align[real_start:real_end][it][1]
				it=t
				while it<=real_end-real_start-1  and conl_align[real_start:real_end][it][0]=="-":
					it=it+1
				if it>real_end-real_start-1:
					a=8
				else:
					a=conl_align[real_start:real_end][it][1]
				print(a,b)
				if int(b)>7 and int(a)>7:
					fuc="no"
			elif w!=v and w!='-' and v=='-':
				delete+=1
				if int(sr[1])>7:
					fuc="no"
				if t==0:
					adj_left=t
				else:
					adj_left=t-1
				bw=wt_seq[real_start:real_end][adj_left]
				bm=record.seq[real_start:real_end][adj_left]
				while adj_left>0 and (bw!=bm or bw==bm=="-"):
					adj_left=adj_left-1
					bw=wt_seq[real_start:real_end][adj_left]
					bm=record.seq[real_start:real_end][adj_left]
				if t==real_end-real_start-1:
					adj_right=t
				else:
					adj_right=t+1
				print(real_start,real_end)
				print(len(wt_seq[real_start:real_end]),t,adj_right)
				aw=wt_seq[real_start:real_end][adj_right]
				am=record.seq[real_start:real_end][adj_right]
				while adj_right<real_end-real_start-1 and (aw!=am or aw==am=="-"):
					adj_right=adj_right+1
					print(adj_right)
					aw=wt_seq[real_start:real_end][adj_right]
					am=record.seq[real_start:real_end][adj_right]
				print(real_start,real_end,adj_left,adj_right)
				if int(conl_align[real_start:real_end][adj_left][1])==9 or int(conl_align[real_start:real_end][adj_right][1])==9:
					fuc="no"
					print(conl_align[real_start:real_end][adj_left])
					print(conl_align[real_start:real_end][adj_right])
			t+=1
		if (con_ch+uncon_ch)>3:
			fuc="no"
		rd[record.id]=(delete,insert,substitute,fuc,str(con_num)+"_"+str(con_ch),str(uncon_num)+"_"+str(uncon_ch),str(record.seq).replace("-",""))
	Fr=open(prefix+"_alleles_deal.xls","w")
	Fr.write("name\t"+"\t".join(df_modified.columns)+"\tp_deleted\tp_inserted\tp_mutated\tp_function\tp_conservation\tp_unconservation\taa_seq\n")
	for index,row in df_modified.iterrows():
		seq_name="temp_%s"%(index)
		if seq_name in rd:
			fm=seq_name+"\t"+"\t".join([str(i) for i in row])+"\t"+"\t".join([str(s) for s in rd[seq_name]])
		else:
			fm=seq_name+"\t"+"\t".join([str(i) for i in row])+"\t-\t-\t-\t-\t-\t-\t-"
		Fr.write(fm+"\n")
	Fr.close()
	#return wt_seq[real_start:real_end]
def my_cut(grna,ref):
	grna=grna.upper()
	ref=ref.upper()
	m=re.search(grna,ref)
	if m:
		cut=m.end()-4 #?
		orientation=1
	else:
		r_grna=str(Seq(grna).reverse_complement())
		m=re.search(r_grna,ref)
		print(r_grna)
		print(ref)
		orientation=-1
		cut=m.start()+3
	return cut,orientation
def get_grna(grna_seq,sample):
	if sample=="20-001" or sample=="20-051":
		grna_seq=grna_seq[8:]
	elif sample=="28-001" or sample=="83-001" or sample=="92-001" or sample=="28-051" or sample=="83-051" or sample=="92-051":
		grna_seq=grna_seq[7:]
	elif sample=="86-001" or sample=="86-051":
		grna_seq=grna_seq[:-4]
	elif sample=="91-001" or sample=="91-051":
		grna_seq=grna_seq[4:]
	return grna_seq
def deal_blast(name,rep,rsp,wt_cds_seq):
	d_align=collections.defaultdict(list)
	d_cds=collections.defaultdict(list)
	for x in open("%s_blast.txt"%name):
		x=x.rstrip()
		l=x.split("\t")
		d_align[l[0]]+=[int(l[2]),int(l[3])]
		d_cds[l[0]]+=[int(l[6]),int(l[7])]
	Fr=open("%s_deal.fa"%name,"w")
	#Fr.write(">wt_cds\n")
	#Fr.write(wt_cds+"\n")
	#Fr.write(">mutant_cds\n")
	#Fr.write(mutant_cds+"\n")
	for record in SeqIO.parse("%s_temp.fa"%name,"fasta"):
		pa=sorted(d_align[record.id])
		pr=sorted(d_cds[record.id])
		if pa:
			start=pa[0]-(pr[0]-int(rsp))
			end=pa[-1]+(int(rep)-pr[-1])
			if start<1:
				start=1
			seq=record.seq[start-1:end]
			Fr.write(">"+record.id+"\n")
			Fr.write(str(seq)+"\n")
		else:
			Fr_align=open("%s_realign.fa"%name,"w")
			Fr_align.write(">wt_cds\n")
			Fr_align.write(wt_cds_seq+"\n")
			Fr_align.write(">"+record.id+"\n")
			Fr_align.write(str(record.seq)+"\n")
			Fr_align.close()
			os.system("mafft --quiet --genafpair %s_realign.fa >%s_realign_mafft.fa"%(name,name))
			mafft_handle=SeqIO.index("%s_realign_mafft.fa"%name,"fasta")
			p=[]
			for i,n in enumerate(mafft_handle["wt_cds"].seq):
				if n!="-":
					p.append(i)
			start=p[0]
			end=p[-1]+1
			seq=str(mafft_handle[record.id].seq[start:end])
			Fr.write(">"+record.id+"\n")
			Fr.write(seq.upper()+"\n")
			#print ">wt_cds"
			#print wt_cds_seq
			#print ">"+record.id
			#print seq
			#sys.exit()
	Fr.close()
def judge_mut(name):
	mafft_handle=SeqIO.index("%s_single_temp_mafft.fa"%name,"fasta")
	wt_seq=str(mafft_handle["wt_cds"].seq)
	mut_seq=str(mafft_handle["mutant_cds"].seq)
	pm=[]
	t=0
	for w,m in zip(wt_seq,mut_seq):
		if w!=m:
			pm.append(t)
			if w!="-" and m!="-":
				return False
		t+=1
	if len(pm)>1:
		return False
	return True

if __name__ == '__main__':
	finish_l=[]
	for ri in glob.glob("*_alleles_deal.xls"):
		finish_l.append(ri.split("_alleles_deal.xls")[0])
	for x in open("first.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		sl=os.popen("find -L . -name \"CRISPResso_on_%s\""%(l[0])).readlines()
		print(sl)
		for p in sl:
			p=p.rstrip()
			name=p.split("/")[1].split("-")[0]+"_"+l[0]
			#if name!="Rep1_48-051":
			#if name in finish_l:
			#	continue
			wt_cds_seq=l[2].upper()
			mutant_cds_seq=l[3].upper()
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
			rep=os.popen("awk '{print $8}' %s_blast.txt|sort|uniq -c|sort -n -r"%name).readlines()[0].rstrip().split()[1]
			rl=os.popen("awk '{print $6}' %s_blast.txt|sort -u"%name).readlines()[0].rstrip()
			if int(rep)!=int(rl):
				print("algn error please check %s_blast.txt"%name)
				sys.exit()
			rsp=os.popen("awk -F'\t' '{print $7}' %s_blast.txt|sort|uniq -c|sort -k 1 -n -r"%name).readlines()[0].rstrip().split()[1]
			print("ref start align position: %s"%rsp)
			#deal_blast(name,rep,rsp,wt_cds_seq)
			#os.system("mafft --globalpair --thread 60 --inputorder %s_deal.fa>%s_temp_mafft.fa"%(name,name))
			grna_seq=get_grna(l[5].upper(),l[0])
			cut_point,orientation=my_cut(grna_seq,wt_cds_seq)
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
				if not judge_mut(name):
					os.system("mafft --quiet --genafpair --thread 3 --inputorder %s_single_temp.fa>%s_single_temp_mafft.fa"%(name,name))
				if not judge_mut(name):
					os.system("mafft --op 0 --ep 1 --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
				if not judge_mut(name):
					os.system("mafft --op 0 --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
				if not judge_mut(name):
					os.system("mafft --globalpair --quiet %s_single_temp.fa >%s_single_temp_mafft.fa"%(name,name))
				mafft_handle=SeqIO.index("%s_single_temp_mafft.fa"%(name),"fasta")
				reald,mut_flag=mut_reads(cut_point,mafft_handle)
				for rni in reald:
					print(rni+"\t"+"".join(reald[rni]))
				print(mafft_handle[record.id].seq)
				#sys.exit()
				#if record.id=="temp_817":
				#	print(reald)
				#	print(mut_flag)
				#	sys.exit()
				translate_and_align(reald,inframed,mut_flag,name,Fr_aa,Fr_log)
			sys.exit()
			Fr_aa.close()
			Fr_log.close()
			faa_num=int(os.popen("wc -l %s.faa"%name).readlines()[0].rstrip().split()[0])
			if faa_num==2:
				os.system("cp %s.faa %s_align.faa"%(name,name))
			else:
				os.system("mafft --quiet --thread 40 --globalpair %s.faa >%s_align.faa"%(name,name))
			aa_file="mut_cds.fa"
			target_aa=l[-1].upper()
			extend=20
			grna_prot=get_grna_prot(wt_cds_seq,grna_seq,extend,cut_point)
			consurf_file=get_consurf_file(target_aa,aa_file)
			#print target_aa
			if not os.path.exists(consurf_file):
				print(consurf_file)
				continue
			get_var(df2,"%s_align.faa"%name,grna_prot,consurf_file,name)
			#break
		#break
