#!/usr/bin/env python
import os,pandas,re,sys,collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import pairwise2
import get_args_and_test
from CRISPResso2 import CRISPRessoShared
def get_cut_point(pa,wt_cds,orientation,grna): #path
	args=get_args_and_test.get_args()
	args.orientation=orientation #need change
	lines=open(pa+"/CRISPResso_RUNNING_LOG.txt").readlines()
	parl=["needleman_wunsch_gap_incentive","quantification_window_center","quantification_window_size","exclude_bp_from_left","exclude_bp_from_right","plot_window_size","needleman_wunsch_aln_matrix_loc","needleman_wunsch_gap_open","needleman_wunsch_gap_extend","aln_seed_len","aln_seed_count","aln_seed_min","default_min_aln_score"]
	for ar in parl:
		m=re.search("--"+ar+" (.*?) ",lines[2].rstrip())
		if not m:
			m=re.search("--"+ar+" (.*?)$",lines[2].rstrip())
		#print ar,m.group(1)
		#print getattr(args,ar)
		if m:				
			if ar!="needleman_wunsch_aln_matrix_loc":
				setattr(args,ar,int(m.group(1)))
			else:
				setattr(args,ar,m.group(1))
		else:
			print ar
			print "===search error===="
			#sys.exit()
		#print getattr(args,ar)
	#sys.exit()
	#m=re.search("-a (\w+?) ",lines[2].rstrip())
	#ref=m.group(1)
	m=re.search("-g (\w+?) ",lines[2].rstrip())
	#grna=m.group(1)
	#print ref
	#print grna
	this_quant_window_coordinates=None
	this_seq=wt_cds
	guides=[grna]
	print getattr(args,"exclude_bp_from_right")
	print getattr(args,"exclude_bp_from_left")
	(this_sgRNA_sequences, this_sgRNA_intervals, this_cut_points, this_sgRNA_plot_offsets, this_include_idxs,this_exclude_idxs, this_plot_idxs) = CRISPRessoShared.get_amplicon_info_for_guides(this_seq,guides,args.quantification_window_center,args.quantification_window_size,this_quant_window_coordinates,args.exclude_bp_from_left,args.exclude_bp_from_right,args.plot_window_size)
	#print this_include_idxs
	print "hehe"
	print args.orientation
	print this_cut_points
	if args.orientation==1:
		cut=this_cut_points[0]
	else:
		cut=this_cut_points[0]+1
	print this_cut_points
	return cut
def get_grna_prot(wt_cds,grna_seq,extend,pa,orientation,Fr_log):
	win=3
	cdsd={}
	for start,end in [(i,i+win) for i in range(0,len(wt_cds),win)]:
		#print start,end-1
		for ti in range(start,end):
			cdsd[ti]=(start,end)
	#print cdsd
	#alignments = pairwise2.align.localxs(wt_cds.upper(),grna_seq.upper(),-1,-1)
	cut=get_cut_point(pa,wt_cds,orientation,grna_seq)
	Fr_log.write("cut_point:%s\n"%cut)
	#start,end=alignments[0][3:]
	start,end=cut-extend,cut+extend
	seq_deal=Seq(wt_cds[cdsd[start][0]:cdsd[end][1]])
	grna_prot=seq_deal.translate()
	#print start,end
	Fr_log.write("cut_extend_region:%s,%s\n"%(start,end))
	#print cdsd[start][0],cdsd[end-1][1]
	Fr_log.write("wt_cds_target_region:%s,%s\n"%(cdsd[start][0],cdsd[end-1][1]))
	#print seq_deal
	Fr_log.write("wt_cds_target_region_seq:%s\n"%(str(seq_deal)))
	#print grna_seq
	#print grna_prot
	return str(grna_prot).upper()
def get_consurf(consurf):
	conl=[]
	for x in open(consurf):
		x=x.rstrip()
		l=x.split("\t")
		conl.append(l)
	return conl
def get_class(saa):
	basicl=["R","H","K"]
	acidl=["D","E"]
	neutral=['A','C','F','G','I','L','M','N','P','Q','S','T','V','W','Y']
	b=0
	a=0
def get_var(blast_file,align_fa,wt_name,grna_prot,Fr_log):
	conl=get_consurf("consurf.txt") #need change	
	seqd=SeqIO.index(align_fa,"fasta")
	wt_seq=str(seqd[wt_name].seq).upper()
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
	wt_seq_deal=wt_seq.replace("-","").upper()
	m=re.search(grna_prot,wt_seq_deal)
	if m:
		start,end=m.span() #end is not include
	else:
		print "grna_prot match error"
		sys.exit()
	#print start,end
	Fr_log.write("wt_target_prot_position:%s,%s\n"%(start,end))
	real_start=pseq[start]
	real_end=pseq[end]
	#print wt_seq[real_start:real_end]
	Fr_log.write("wt_target_align_prot_position:%s,%s\n"%(real_start,real_end))
	#sys.exit()
	Fr_log.write("conl_align_seq:%s\n"%"".join([i[0] for i in conl_align[real_start:real_end]]))
	Fr_log.write("conl_align_seq_score:%s\n"%"".join([str(i[1]) for i in conl_align[real_start:real_end]]))
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
		for w,v,sr in zip(wt_seq[real_start:real_end],record.seq[real_start:real_end],conl_align[real_start:real_end]):
			saa=sr[2].split(",")
			#clas=get_class(saa)
			if int(sr[1])>5:
				con_num+=1
			else:
				uncon_num+=1
			if w!=v:
				if int(sr[1])>5:
					con_ch+=1
				else:
					uncon_ch+=1
			if v=="*":
				fuc="no"
			if w!=v and w!="-" and v!="-":
				substitute+=1
				if int(sr[1])>5 and (not v in saa):
					fuc="no"
				#else:
					#if (not v in saa) and (not v in clas):
				#	if not v in saa:
				#		fuc="no"
			elif w!=v and w=='-' and v!="-":
				insert+=1
				print conl_align[real_start:real_end][t-1][1]
				print conl_align[real_start:real_end][t+1][1]
				if int(conl_align[real_start:real_end][t-1][1])>5 and int(conl_align[real_start:real_end][t+1][1])>5:
					fuc="no"
			elif w!=v and w!='-' and v=='-':
				delete+=1
				if int(sr[1])>5:
					fuc="no"
			t+=1
		rd[record.id]=(delete,insert,substitute,fuc,str(con_num)+"_"+str(con_ch),str(uncon_num)+"_"+str(uncon_ch))
	lines=open(blast_file).readlines()
	Fr=open(blast_file.split(".")[0]+"_deal.xls","w")
	Fr.write(lines[0].rstrip()+"\tp_deleted\tp_inserted\tp_mutated\tp_function\tp_conservation\tp_unconservation\n")
	for x in lines[1:]:
		x=x.rstrip()
		l=x.split("\t")
		fm="\t".join([str(s) for s in rd[l[0]]])
		Fr.write(x+"\t"+fm+"\n")
	Fr.close()
	return wt_seq[real_start:real_end]
def get_matrix(align_fa,name,blast_file):
	dc={}
	for x in open(blast_file):
		x=x.rstrip()
		l=x.split("\t")
		dc[l[0]]=l[-2]
	aa=list(IUPAC.protein.letters)+["-"]
	d=collections.defaultdict(dict)
	for record in SeqIO.parse(align_fa,"fasta"):
		if record.id=="wt":
			wt_aa=str(record.seq)
			continue
		for i in range(len(record.seq)): 
			d[i][record.seq[i]]=d[i].get(record.seq[i],0)+int(dc[record.id])
	#print d
	conl=get_consurf("consurf.txt")
	conl_align=[]
	t=0
	for n,s in enumerate(wt_aa):
		if s!="-":
			conl_align.append(conl[t])
			t+=1
		else:
			conl_align.append(["-",0,"0"])
	Fr=open(name+"_aa_matrix.xls","w")
	Fr.write("aa\t"+"\t".join([str(i) for i in range(len(d.keys()))])+"\n")
	Fr.write("wt\t"+"\t".join(list(wt_aa))+"\n")
	Fr.write("conser\t"+"\t".join([str(i[1]) for i in conl_align])+"\n")
	for a in aa:
		fm=a+"\t"
		for i in range(len(d.keys())):
			if a in d[i]:
				fm+=str(d[i][a])+"\t"
			else:
				fm+="0\t"
		Fr.write(fm[:-1]+"\n")
	Fr.close()
def get_final_results(pa_file,wt_cds,wt_aa,grna_seq,extend,orientation,Fr_log):
	pa=os.path.dirname(pa_file)
	name=pa.split("/")[1].replace("CRISPResso_on_","")
	print name
	Fr_log.write("sample:%s\n"%name)
	grna_prot=get_grna_prot(wt_cds,grna_seq,extend,pa,orientation,Fr_log)
	Fr=open("%s.fa"%name,"w")
	#df=pandas.read_csv("Alleles_frequency_table.txt",sep="\t",header=0,index_col=None)
	df=pandas.read_csv(pa_file,sep="\t",header=0,index_col=None) #need change
	df2=df[(df["Read_Status"]=="MODIFIED")&((df["n_inserted"]-df["n_deleted"]+1)%3==0)] #need change
	d={}
	for index,row in df2.iterrows():
		seq=row["Aligned_Sequence"]
		seq=seq.replace("-","")
		Fr.write(">temp_%s\n"%(index))
		Fr.write(seq+"\n")
		d["temp_%s"%index]=row["Aligned_Sequence"]+"\t"+row["Reference_Sequence"]+"\t"+row["Reference_Name"]+"\t"+row["Read_Status"]+"\t"+str(row["n_deleted"])+"\t"+str(row["n_inserted"])+"\t"+str(row["n_mutated"])+"\t"+str(row["#Reads"])+"\t"+str(row["%Reads"])
	Fr.close()
	os.system("blastn -db cds -query %s.fa -out %s_blast.out -num_threads 40 -evalue 1e-10 -outfmt \"6 qseqid qlen qstart qend sseqid slen sstart send qcovs bitscore evalue pident\""%(name,name))
	Fr=open("%s_blast.xls"%name,"w")
	Fr.write("query\tqlen\tqstart\tqend\ttarget\ttlen\ttstart\ttend\tqcov\tscore\tevalue\tidentify\tAligned_Sequence\tReference_Sequence\tReference_Name\tRead_Status\tdeleted\tinserted\tmutated\t#Reads\t%Reads\n")
	for x in open("%s_blast.out"%name):
		x=x.rstrip()
		l=x.split("\t")
		Fr.write(x+"\t"+d[l[0]]+"\n")
	Fr.close()
	#al=os.popen("awk '{print $8-$7+1}' %s_blast.out|sort|uniq -c|sort -n -r"%name).readlines()[0].rstrip().split()[1]
	rep=os.popen("awk '{print $8}' %s_blast.out|sort|uniq -c|sort -n -r"%name).readlines()[0].rstrip().split()[1]
	rl=os.popen("awk '{print $6}' %s_blast.out|sort -u"%name).readlines()[0].rstrip()
	if int(rep)!=int(rl):
		print "algn error please check %s_blast.out"%name
		sys.exit()
	rsp=os.popen("awk -F'\t' '{print $7}' %s_blast.out|sort|uniq -c|sort -k 1 -n -r"%name).readlines()[0].rstrip().split()[1]
	print "ref start align position: %s"%rsp
	d_align=collections.defaultdict(list)
	d_cds=collections.defaultdict(list)
	for x in open("%s_blast.out"%name):
		x=x.rstrip()
		l=x.split("\t")
		d_align[l[0]]+=[int(l[2]),int(l[3])]
		d_cds[l[0]]+=[int(l[6]),int(l[7])]
	Fr=open("%s_deal.fa"%name,"w")
	Fr_aa=open("%s_deal.faa"%name,"w")
	Fr_aa.write(">wt\n") #wt_name
	Fr_aa.write(wt_aa+"\n")
	for record in SeqIO.parse("%s.fa"%name,"fasta"):
		pa=sorted(d_align[record.id])
		pr=sorted(d_cds[record.id])
		abn=0
		if len(pa)!=2:
			abn=1
		if pr[0]!=int(rsp):
			abn=1
			start=pa[0]-(pr[0]-int(rsp))
			#end=pa[-1]
			if start<1:
				start=1
		else:
			start=pa[0]
		if pr[-1]!=int(rep):
			abn=1
			end=pa[-1]+(int(rep)-pr[-1]) #may be need change
		else:
			end=pa[-1]
		seq=record.seq[start-1:end]
		if abn!=0:
			Fr.write(">"+record.id+" abnormal\n")
			Fr_aa.write(">"+record.id+" abnormal\n")
		else:
			Fr.write(">"+record.id+"\n")
			Fr_aa.write(">"+record.id+"\n")
		Fr.write(str(seq)+"\n")
		Fr_aa.write(str(seq.translate())+"\n")
		#break
	Fr.close()
	Fr_aa.close()
	os.system("mafft --auto --reorder %s_deal.faa >%s_deal_align.faa"%(name,name))
	get_matrix("%s_deal_align.faa"%name,name,"%s_blast.xls"%name)
	grna_prot_align=get_var("%s_blast.xls"%name,"%s_deal_align.faa"%name,"wt",grna_prot,Fr_log)
	Fr_log.write("wt_cds:%s\n"%wt_cds)
	Fr_log.write("wt_aa:%s\n"%wt_aa)
	Fr_log.write("grna_seq:%s\n"%grna_seq)
	Fr_log.write("target_region_grna_prot:%s\n"%grna_prot)
	Fr_log.write("grna_prot_align:%s\n"%grna_prot_align)
	get_func_ratio("%s_blast_deal.xls"%name,Fr_log)
def get_func_ratio(blast_deal,Fr_log):
	df=pandas.read_csv(blast_deal,sep="\t",header=0,index_col=0)
	yes=df[df["p_function"]=='yes']["#Reads"].sum()
	al=df["#Reads"].sum()
	Fr_log.write("%s yes_ration:%s\n"%(blast_deal,yes/float(al)))
if __name__ == '__main__':
	lines=os.popen("find . -name 'Pure_Alleles_frequency_table.txt'").readlines() #need change
	wt_cds="GACGGCAAACTGCTCGATATCAATAAAGACTTCCAGCCGTATTACGGGGAAGGAGGGCGCATTCTGGAGATTCGGACACCTGAGGCAGTGACGAGCATCAAGAAGCGAGGAGAAAGCTTGGGGTACACAGAAGGGGCCTTGCTGGCCTTGGCCTTCATCATCATCCTCTGTTGCATCCCAGCCATCTTGGTCGTCTTAGTAAGCTACCGACA"
	wt_aa="DGKLLDINKDFQPYYGEGGRILEIRTPEAVTSIKKRGESLGYTEGALLALAFIIILCCIPAILVVLVSYR"
	grna_seq='AGCGAGGAGAAAGCTT'
	orientation=1
	Fr_log=open("aa.log","w")
	extend=20
	for x in lines:
		x=x.rstrip()
		print x
		#os.system("unzip -o %s"%x) 
		#l=x.split("/")
		get_final_results(x,wt_cds,wt_aa,grna_seq,extend,orientation,Fr_log)
		#break
	Fr_log.close()
