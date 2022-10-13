import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


def StatsInsertSizes(InsertList, outpath, title="Size of inserts accross Samples", sampleList=None,fontsize=80,outtxt=None, threshold=(None,None), outliers=None ):
	"""StatsInsertSizes(InsertList, outpath, title="Size of inserts accross Samples", sampleList=None)---->
	InsertList: list of files containing Insert metrics, outputs of picard.analysis.InsertSizeMetrics
	outpath: path to output image. Extension determines format (svg,png,jpg,...)
	title: Title to give to figure
	sampleList: list of sample names or alias. If given, should be part of insert file's names.
	outtxt: path to output results as TSV file
	threshold: couple of numeric values (tmin, tmax). If outliers is set, will write to outliers sample names that have mean insert size below tmin or above tmax
	outliers: path to file where outlier's names should be written

	This function synthetizes insert size metrics performed by picard in a barplot representing the mean insert size as well as its standard deviation in each sample"""
	Lab=InsertList
	Results=[]
	if sampleList!=None:
		Lab=sampleList
	for i in range(len(Lab)):
		l=Lab[i]
		print("Processing "+l+" ("+str(i+1)+"/"+str(len(Lab))+")" )
		file=None
		for f in InsertList:
			if l in f: 
				file=f
				break
		stats=open(file,'r')
		foundmean=False
		while not foundmean:
			line=stats.readline()
			if "METRICS CLASS" in line:
				foundmean=True
		t=(stats.readline()).split('\t')
		v=(stats.readline()).split('\t')
		meanIdx=t.index("MEAN_INSERT_SIZE")
		stdIdx=t.index("STANDARD_DEVIATION")
		Results.append([l,v[meanIdx],v[stdIdx]])
	if outtxt!=None:
		with open(outtxt,"w") as txt:
			txt.write("\n".join((["Samples\tMean\tStd"]+["\t".join(l) for l in Results])))
	Results=np.array(Results)
	#Results=np.array(Results, dtype=[("label","str"),("mean","float"),("std","float")])
	Results=Results[Results[:,1].astype('float').argsort()]
	if outliers!=None:
		with open(outliers,'a') as outliersFile:
			tmin= threshold[0] if threshold[0]!=None else 0
			tmax= threshold[1] if threshold[1]!=None else max(Results[:,1])
			Oultiers_min= Results[Results[:,1]<tmin,0]
			Outliers_max= Results[Results[:,1]>tmax,0]
			for om in Oultiers_min:
				outliersFile.write(om+'\t#Inserts too short')
			for oM in Outliers_max:
				outliersFile.write(oM+'\t#Inserts too long')
	fig,ax=plt.subplots(figsize=(250,40))
	plt.xticks(rotation=90)
	ax.set_xticks(range(len(Lab)))
	ax.set_xticklabels(Lab)
	ax.bar(range(len(Lab)), Results[:,1].astype('float'), yerr=Results[:,2].astype('float'))
	plt.title(title,fontsize=fontsize)
	plt.xlabel("Samples",fontsize=fontsize)
	plt.ylabel("Mean insert size",fontsize=fontsize)
	plt.savefig(outpath)