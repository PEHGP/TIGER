#!/usr/bin/env python2.7
import sys,pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
f=sys.argv[1]
df=pandas.read_csv(f,sep="\t",header=0,index_col=[0,1,9])
aa="aaxd3"
col=[aa,"%Inframe","%Indel","Specificity_Score","%Function_in_inframe"]
df=df.loc[:,col]
df.sort_values(by=['sort_grna'],inplace=True)
theta = np.linspace(0, 360,len(df.columns), endpoint=False)
adj_angle = theta[-1] + 90 - 360
theta += adj_angle
#print theta
X_ticks = np.radians(theta)
X = np.append(X_ticks,X_ticks[0])
for tittle,gdf in df.groupby(df.index.get_level_values(0)):
	num=9 #need change,Affect image size,If the number of grna is greater than num, you need to change the value of num 
	fig, axs = plt.subplots(1,num,figsize=(60,20),dpi=300,subplot_kw=dict(projection='polar'))
	ym=100
	j=0
	print tittle
	lab=[]
	gdf2=gdf.copy()
	jj=gdf2.shape[0]
	for index,row in gdf2.iterrows():
		print index
		row["%Indel"]=row["%Indel"]/10*ym
		row["%Inframe"]=row["%Inframe"]/100*ym
		#row["%No AA Indel"]=row["%No AA Indel"]/100*ym
		#row["Distance to Mutation"]=(15-row["Distance to Mutation"])/15*ym
		row[aa]=float(row[aa])/100*ym
		row["%Function_in_inframe"]=float(row["%Function_in_inframe"])/100*ym
		row["Specificity_Score"]=row["Specificity_Score"]/55*ym
		Y=np.array(row)
		Y=np.append(Y,Y[0])
		axs[j].plot(X,Y,marker='')
		axs[j].fill(X,Y,alpha=0.25)
		#lab.append(index)
		axs[j].set_ylim(0,ym)
		axs[j].set_xticks(X)
		axs[j].set_xticklabels(gdf2.columns)
		axs[j].set_yticklabels([])
		axs[j].tick_params(rotation='auto')
		axs[j].grid(True)
		axs[j].set_title(index[1])
		if jj!=1:		
			axs[jj].plot(X,Y,marker='')
			axs[jj].fill(X,Y,alpha=0.25)
			lab.append(index[1])
		j+=1
	if jj!=1:
		axs[jj].set_ylim(0,ym)
		axs[jj].set_xticks(X)
		axs[jj].set_xticklabels(gdf2.columns)
		axs[jj].set_yticklabels([])
		axs[j].tick_params(rotation='auto')
		axs[jj].grid(True)
		axs[jj].set_title(tittle)
		axs[jj].legend(lab,bbox_to_anchor=(1.17, 1),loc='upper right')
	gdf3=gdf.copy()
	if gdf3.shape[0]!=1:
		#gdf3["Distance to Mutation"]=gdf3["Distance to Mutation"].max()-gdf3["Distance to Mutation"]
		#gdf3=gdf3/gdf3.max()*ym
		#gdf3=(gdf3-gdf3.min())/(gdf3.max()-gdf3.min())*ym
		#gdf3.fillna(0,inplace=True)
		#print gdf
		lab=[]
		j+=1
		for index,row in gdf3.iterrows():
			print index
			if row[aa]=="-":
				continue
			Y=np.array(row)
			Y=np.append(Y,Y[0])
			axs[j].plot(X,Y,marker='')
			axs[j].fill(X,Y,alpha=0.25)
			lab.append(index[1])
		axs[j].set_ylim(0,ym)
		axs[j].set_xticks(X)
		axs[j].set_xticklabels(gdf3.columns)
		axs[j].set_yticklabels([])
		axs[j].tick_params(rotation='auto')
		axs[j].grid(True)
		axs[j].set_title(tittle)
		axs[j].legend(lab,bbox_to_anchor=(1.17, 1),loc='upper right')
	if jj!=1:
		j+=1
	axs[j].set_ylim(0,ym)
	axs[j].set_xticks(X)
	axs[j].set_xticklabels(gdf3.columns)
	#axs[j].set_yticklabels([])
	axs[j].tick_params(rotation='auto')
	axs[j].grid(True)
	axs[j].set_title("legend")
	fig.tight_layout()
	fig.savefig("%s_polar.svg"%tittle,format='svg')
	fig.clf()
	plt.close(fig)
	#break
