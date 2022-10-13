import os
import sys
from snakemake.utils import *


WF_dir=workflow.basedir
config_dir=os.path.dirname(WF_dir)
sys.path.insert(1, os.path.join(WF_dir,"scripts"))
from Utils import StatsInsertSizes


configfile: os.path.join(config_dir,"config/config.yaml")

#config_dir=os.path.dirname(os.path.dirname(str(workflow.configfile)))

validate(config,os.path.join(config_dir,"config/config.schema.yaml"))


prefix=config["prefix"]
CRAM_list=config["CRAM_list"]
samples_list=config["samples_list"]
tempdir=config["tempdir"]
outdir=config["outdir"]

insertcompute_env=os.path.join(WF_dir,'envs/insertcompute.yaml')
statsinsert_env=os.path.join(WF_dir,'envs/statsinsert.yaml')
samrun_env=os.path.join(WF_dir,'envs/samrun.yaml')
somalier_env=os.path.join(WF_dir,'envs/somalier.yaml')


workdir: tempdir

SAMPLES=[]
with open(samples_list,'r') as samp_file:
	for s in samp_file.readlines():
		SAMPLES.append(s.replace('\n',''))

CRAMS=[]
with open(CRAM_list,'r') as cram_file:
	for s in cram_file.readlines():
		CRAMS.append(s.replace('\n',''))

# module coverage:
# 	snakefile:
# 		os.path.join(WF_dir,"rules/GetCover.smk")
# 	config:
# 		config

# use rule * from coverage as cov_*



rule all:
	input:
		txt=os.path.join(outdir,"Insert/"+prefix+".insert_size_metrics.tsv"),
		InOff=os.path.join(outdir,"Stats/"+prefix+".INOFF.tsv"),
		stats=os.path.join(outdir,"Stats/"+prefix+".Stats.tsv"),
		#tsv=rules.cov_all.input.tsv,
		relate=os.path.join(outdir,"Relatedness/"+prefix+".pairs.tsv")


rule InsertCompute:
	input: 
		cram=CRAMS
	output: 
		pdf=temp("{sample}.insert_size_histogram.pdf"),
		txt=temp("{sample}.insert_size_metrics.txt"),
	params:
		reference=config["reference"],
	threads:
		1
	conda:
		insertcompute_env
	shell:
		"""
		SAMPLE=`samtools samples {input.cram}|cut -f1`
		picard CollectInsertSizeMetrics -I {input.cram} -O ${{SAMPLE}}.insert_size_metrics.txt -H ${{SAMPLE}}.insert_size_histogram.pdf -M 0.5 -R {params.reference}
		"""

rule StatsInserts:
	input:
		txt=expand("{sample}.insert_size_metrics.txt",sample=SAMPLES)
	output:
		svg=temp(os.path.join(outdir,"Insert/"+prefix+".insert_size_metrics.svg")),
		txt=os.path.join(outdir,"Insert/"+prefix+".insert_size_metrics.tsv"),
	params:
		samples=SAMPLES
	run:
		StatsInsertSizes(input.txt, output.svg, title="Size of inserts accross Samples", sampleList=params.samples,outtxt=output.txt)


rule SamRun:
	input:
		cram=CRAMS,
	output:
		InOff=temp("{sample}.INOFF.tsv"),
		Stats=temp("{sample}.Stats.tsv")
	params:
		reference=config["reference"],
		regions=config["regions"],
		mapq=30
	conda:
		samrun_env
	shell:
		"""
		SAMPLE=`samtools samples {input.cram} | cut -f1`
		TOT=`samtools view -c -q {params.mapq} --reference {params.reference} {input.cram}`
		IN=`samtools view -c -q {params.mapq} -L {params.regions} --reference {params.reference} {input.cram}`
		echo -e '${{IN}}\t${{TOT}}'>{output.InOff}

		samtools flagstats -O tsv {input.cram} > {output.Stats}
		"""

rule catStats:
	input:
		InOff=expand("{sample}.INOFF.tsv",sample=SAMPLES),
		Stats=expand("{sample}.Stats.tsv",sample=SAMPLES)
	output:
		InOff=os.path.join(outdir,"Stats/"+prefix+".INOFF.tsv"),
		Stats=os.path.join(outdir,"Stats/"+prefix+".Stats.tsv"),
	run:
		INOFF=[['sample','IN-target','OFF-target']]
		STATS=[]
		for inOff in input.InOff:
			samp=inOff.split('.')[-3]
			with open(inOff) as inOffFile:
				for l in inOffFile.readlines():
					INOFF.append([samp]+l.replace('\n','').split('\t'))

		for stats in input.Stats:
			samp=stats.split('.')[-3]
			with open(stats) as statsFile:
				s=[samp]
				titles=['sample']
				for l in statsFile.readlines():
					L=l.replace('\n','').split('\t')
					s.append(L[0])
					titles.append(L[-1])
				STATS.append(s)
		STATS=[titles]+STATS

		with open(output.InOff,'w') as OUT_InOff:
			OUT_InOff.write('\n'.join(['\t'.join(l) for l in INOFF]))

		with open(output.Stats,'w') as OUT_Stats:
			OUT_Stats.write('\n'.join(['\t'.join(l) for l in STATS]))



rule somalier_extract:
	input:
		cram=CRAMS
	output:
		som=temp("{sample}.somalier")
	conda:
		somalier_env
	params:
		reference=config["reference"],
		sites=config["sites_som"]
	shell:
		"somalier extract --sites {params.sites} -f {params.reference} {input.cram}"



rule somalier_relate:
	input:
		som=expand("{sample}.somalier",sample=SAMPLES)
	output:
		relate=os.path.join(outdir,"Relatedness/"+prefix+".pairs.tsv")
	conda:
		somalier_env
	params:
		outprefix=os.path.join(outdir,"Relatedness/"+prefix)
	shell:
		"somalier relate  -o {params.outprefix} {input.som}"
