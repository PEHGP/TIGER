#!/usr/bin/env python
import sys,re
import inDelphi
from Bio.Seq import Seq
def my_cut(grna,ref):
	grna=grna.upper()
	ref=ref.upper()
	m=re.search(grna,ref)
	if m:
		cut=m.end()-3 #?
		orientation=1
	else:
		ref=str(Seq(ref).reverse_complement())
		m=re.search(grna,ref)
		orientation=-1
		cut=m.end()-3
	return ref,cut,orientation
if __name__ == '__main__':
	f=sys.argv[1] #input
	inDelphi.init_model(celltype='U2OS')
	lines=open(f).readlines()
	h=lines[0].rstrip()+"\twt_cut\twt_orientation\tmut_cut\tmut_orientation"
	Fr=open(f.split(".")[0]+"_ori.xls","w")
	Fr.write(h+"\n")
	for x in lines[1:]:
		x=x.rstrip()
		l=x.split("\t")
		wt_grna=l[2].upper()
		mut_grna=l[3].upper()
		wt_ref=l[6].upper()
		mut_ref=l[7].upper()
		wt_ref,wt_cut,wt_orientation=my_cut(wt_grna,wt_ref)
		mut_ref,mut_cut,mut_orientation=my_cut(mut_grna,mut_ref)
		#print wt_cut,wt_orientation
		#print mut_cut,mut_orientation
		fm="\t".join(l[:5])+"\t"+wt_ref+"\t"+mut_ref+"\t"+str(wt_cut)+"\t"+str(wt_orientation)+"\t"+str(mut_cut)+"\t"+str(mut_orientation)
		Fr.write(fm+"\n")
		wt_pred_df, wt_stats = inDelphi.predict(wt_ref,wt_cut)
		mut_pred_df, mut_stats = inDelphi.predict(mut_ref,mut_cut)
		#wt_stats['gRNA']=wt_grna
		#mut_stats['gRNA']=mut_grna
		#wt_stats['gRNA orientation']=wt_orientation
		#mut_stats['gRNA orientation']=mut_orientation
		#print wt_pred_df
		wt_df_indel=inDelphi.get_indel_length_fqs(wt_pred_df)
		mut_df_indel=inDelphi.get_indel_length_fqs(mut_pred_df)
		wt_df_indel.to_csv("wt_"+l[1]+"_indel.xls",sep="\t",index=False)
		mut_df_indel.to_csv("mut_"+l[1]+"_indel.xls",sep="\t",index=False)
		#break
	Fr.close()
