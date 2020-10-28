#!/usr/bin/env python
import collections,os
d=collections.defaultdict(list)
for x in open("200115-input.txt").readlines()[1:]:
	x=x.rstrip()
	l=x.split("\t")
	d[l[0]].append("\t".join(l[1:]))
for ind in d:
	Fr=open("Index%s_amplicons.txt"%ind,"w")
	for li in d[ind]:
		Fr.write(li+"\n")
	Fr.close()
	cmd="CRISPRessoPooled -q 30 -o Index%s_Output -w 12 --exclude_bp_from_left 10 --exclude_bp_from_right 10 --min_frequency_alleles_around_cut_to_plot 0 --max_rows_alleles_around_cut_to_plot 200 --crispresso1_mode --suppress_report --max_paired_end_reads_overlap 150 -f Index%s_amplicons.txt -n Index%s -r1 Index%s_combined_R1.fastq.gz -r2 Index%s_combined_R2.fastq.gz"%(ind,ind,ind,ind,ind)
	print cmd
	os.system(cmd)
	#break
