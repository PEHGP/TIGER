#!/usr/bin/env python2.7
import os,re,sys
lines=os.popen("find . -name \"*.zip\"").readlines()
stat=["sample\tAlleles_frequency_table\tPure_Alleles_frequency_table\tAlleles_frequency_table_around_sgRNA\tPure_Alleles_frequency_table_around_sgRNA"]
for x in lines:
	x=x.rstrip()
	os.system("unzip -o %s"%x)
	if "WT-genome".upper() in x.upper():
		r="GAAGAAAG"
	elif "3j-genome".upper() in x.upper():
		r="GA-GAAAG"
	else:
		print "=============error============"
		print x
		sys.exit()
	p=os.path.dirname(os.path.abspath(x))
	#print p
	os.system("awk -F'\t' '$1!~/%s/||$1==Aligned_Sequence{print $0}' Alleles_frequency_table.txt >Pure_Alleles_frequency_table.txt"%r)
	os.system("awk -F'\t' '$1!~/%s/||$1==Aligned_Sequence{print $0}' %s/Alleles_frequency_table_around_sgRNA_*.txt >Pure_Alleles_frequency_table_around_sgRNA.txt"%(r,p))
	
	at=os.popen("wc -l Alleles_frequency_table.txt").readlines()[0].rstrip().split()[0]
	pat=os.popen("wc -l Pure_Alleles_frequency_table.txt").readlines()[0].rstrip().split()[0]
	atg=os.popen("wc -l %s/Alleles_frequency_table_around_sgRNA_*.txt"%p).readlines()[0].rstrip().split()[0]
	patg=os.popen("wc -l Pure_Alleles_frequency_table_around_sgRNA.txt").readlines()[0].rstrip().split()[0]
	os.system("mv Pure_Alleles_frequency_table.txt %s/"%(p))
	os.system("mv Pure_Alleles_frequency_table_around_sgRNA.txt %s/"%(p))
	os.system("rm -rf Alleles_frequency_table.txt")
	stat.append(x+"\t"+at+"\t"+pat+"\t"+atg+"\t"+patg)
	#break
Fr=open("stat.xls","w")
Fr.write("\n".join(stat)+"\n")
Fr.close()
