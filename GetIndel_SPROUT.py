#!/usr/bin/env python
import numpy as np
import pickle,re
import sys,pandas
from Bio import SeqIO
from Bio.Seq import Seq
def one_hot_index(nucleotide):
	if nucleotide == 'g':
		nucleotide = 'G'
	elif nucleotide == 'a':
		nucleotide = 'A'
	elif nucleotide == 'c':
		nucleotide = 'C'
	elif nucleotide == 't':
		nucleotide = 'T'
	nucleotide_array = ['A', 'C', 'G', 'T']
	return nucleotide_array.index(nucleotide)
if __name__ == '__main__':
	m_frac_total_ins = pickle.load(open('models/fraction_total_insertions_other_cells.p', 'rb'))
	m_frac_total_del = pickle.load(open('models/fraction_total_deletions_other_cells.p', 'rb'))
	m_frac_mutant_ins = pickle.load(open('models/fraction_insertions_other_cells.p', 'rb'))
	m_avg_ins_length = pickle.load(open('models/exp_ins_length_other_cells.p', 'rb'))
	m_avg_del_length = pickle.load(open('models/exp_deletion_length_other_cells.p', 'rb'))
	m_diversity = pickle.load(open('models/diversity_other_cells.p', 'rb'))
	single_bp_inserted = pickle.load(open('models/single_insertion_type_4class-classification_other_cells.p', 'rb'))
	dr={}
	for x in open("200709-inDelphi-input.txt").readlines()[1:]:
		x=x.rstrip()
		l=x.split("\t")
		grna=l[2].upper()
		ref=l[6].upper()
		name=l[1]
		dr[name]={}
		m=re.search(grna,ref)
		if m:
			start=m.start()
			end=m.end()
			#print(grna)
			#print(ref[start:end])
		else:
			ref=str(Seq(ref).reverse_complement())
			m=re.search(grna,ref)
			start=m.start()
			end=m.end()
			#print(grna)
			#print(ref[start:end])
		grna=ref[start-3:end+3]
		print(len(grna),grna)
		spacer_pam=grna
		sequence_pam_per_gene_grna = np.zeros((1, 23, 4), dtype=bool)
		for ind,basepair in enumerate(spacer_pam):
			sequence_pam_per_gene_grna[0,ind,one_hot_index(basepair)] = 1
		#print(sequence_pam_per_gene_grna)
		sequence_pam_per_gene_grna = np.reshape(sequence_pam_per_gene_grna , (1,-1))
		#print(sequence_pam_per_gene_grna)
		frac_total_ins=float(m_frac_total_ins.predict(sequence_pam_per_gene_grna)[0]) #Fraction of total reads with insertion
		frac_mutant_ins=float(m_frac_mutant_ins.predict(sequence_pam_per_gene_grna)[0]) #fraction of indel mutant reads with an insertion
		Inse_del_ratio=frac_mutant_ins/(1 -frac_mutant_ins) #Insertion to deletion ratio
		avg_ins_length=float(m_avg_ins_length.predict(sequence_pam_per_gene_grna)[0])#Average insertion length
		avg_del_length=float(m_avg_del_length.predict(sequence_pam_per_gene_grna)[0])#Average deletion length
		diversity = m_diversity.predict(sequence_pam_per_gene_grna)[0]
		nucleotide_array = ['A', 'C', 'G', 'T']
		ins_base=nucleotide_array[int(single_bp_inserted.predict(sequence_pam_per_gene_grna)[0])]
		dr[name]["frac_total_ins"]=frac_total_ins
		dr[name]["Inse_del_ratio"]=Inse_del_ratio
		dr[name]["avg_ins_length"]=avg_ins_length
		dr[name]["avg_del_length"]=avg_del_length
		dr[name]["diversity"]=avg_del_length
		dr[name]["ins_base"]=ins_base
		#break
	df=pandas.DataFrame.from_dict(dr,orient='index')
	df.to_csv("SPROUT_result.xls",sep="\t",index_label='grna')

