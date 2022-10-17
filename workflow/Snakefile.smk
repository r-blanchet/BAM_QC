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

test_path=''
is_test=config['test']
if is_test:
	test_path=config_dir
	with open(os.path.join(config_dir,config["CRAM_list"]),'w') as cramlist:
		cramlist.write(os.path.join(config_dir,'.test/data_test/Test.bam'))


prefix=config["prefix"]
CRAM_list=os.path.join(test_path,config["CRAM_list"])
samples_list=os.path.join(test_path,config["samples_list"])
reference=os.path.join(test_path,config['reference'])
sites_som=os.path.join(test_path,config['sites_som'])
regions=os.path.join(test_path,config['regions'])
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



rule all:
	input:
		txt=os.path.join(outdir,"Insert/"+prefix+".insert_size_metrics.tsv"),
		InOff=os.path.join(outdir,"Stats/"+prefix+".INOFF.tsv"),
		stats=os.path.join(outdir,"Stats/"+prefix+".Stats.tsv"),
		relate=os.path.join(outdir,"Relatedness/"+prefix+".pairs.tsv")


rule InsertCompute:
	input: 
		cram=CRAMS
	output: 
		pdf=temp("{sample}.insert_size_histogram.pdf"),
		txt=temp("{sample}.insert_size_metrics.txt"),
	params:
		reference=reference,
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
		reference=reference,
		regions=regions,
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
		reference=reference,
		sites=sites_som
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
