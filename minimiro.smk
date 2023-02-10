import os 
import sys
import re
import re
import pandas as pd 



configfile: "minimiro.yaml"
SDIR=os.path.dirname(workflow.snakefile)

shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")

SMS = list(config.keys())
SEQS=["ref", "query"]
SCORES=config["scores"];  SMS.remove("scores") # minimum peak dp socre 

RS = {}
RGNS = {}
QS = {}
QRGNS = {}
RCS = {}
GENES = {}
BAMS = {}
for SM in SMS:
	RS[SM] = config[SM]["ref"] # reference seqeunce 
	assert os.path.exists(RS[SM]+".fai")
	RGNS[SM]	=	config[SM]["regions"] # region(s) to pull from reference 
	QS[SM]		=	config[SM]["query"] # query sequence
	assert os.path.exists(QS[SM]+".fai")
	QRGNS[SM]	=	config[SM]["queryregions"] # region(s) to pull from query
	if("rc" in config[SM]):
		RCS[SM] = config[SM]["rc"] 
	else: 
		RCS[SM] = False
	GENES[SM] = config[SM]["genes"]	

	if("bam" in config[SM]):
		BAMS[SM] = config[SM]["bam"]

wildcard_constraints:
	SEQ="|".join(SEQS),
	SM="|".join(SMS),
	SCORE="\d+",

rule all:
	input:
		pdf	= expand("minimiro/{SM}_{SCORE}_aln.pdf", SM=SMS, SCORE=SCORES),
		png	= expand("minimiro/{SM}_coverage.png", SM=list(BAMS.keys())),


# -------- Input Functions -------- #

def get_ref(wildcards):
	SM = str(wildcards.SM)
	return(RS[SM])

def get_query(wildcards):
	SM = str(wildcards.SM)
	return(QS[SM])

def get_ref_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( RGNS[SM] ) )

def get_query_rgns(wildcards):
	SM = str(wildcards.SM)
	return( " ".join( QRGNS[SM] ) )

def get_rc(wildcards):
	SM = str(wildcards.SM)
	return(RCS[SM])

def get_genes(wildcards):
	SM = str(wildcards.SM)
	return(GENES[SM])

def get_ref_bed(wildcards):
	SM = str(wildcards.SM)
	rtn = []
	for rgn in RGNS[SM]:
		match = re.match("(.+):(\d+)-(\d+)", rgn.strip())
		if(match):
			rtn.append('{}\t{}\t{}\n'.format(*match.groups()))
		else:
			rtn.append('{}\t{}\t{}\n'.format(rgn, 0, 1000000000))
	return( rtn )

def get_score(wildcards):
	return( int(str(wildcards.SCORE)))

def get_bam(wc):
	SM = str(wc.SM)
	return(BAMS[SM])

def get_rb_msize(wildcards):
	SM = str(wildcards.SM)
	return int(config[SM].get('rb_msize', 1000))


# -------- Begin rules -------- #

rule get_rgns:
	input:
		ref=get_ref, 
		query=get_query,
	output:
		ref = temp("temp/{SM}_ref.fasta"),
		query = temp("temp/{SM}_query.fasta"),
	params:
		rgns = get_ref_rgns,
		qrgns = get_query_rgns,
		rc = get_rc,
	threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (1 * threads),
        hrs = 72
	run:
		shell("samtools faidx {input.ref} {params.rgns} > {output.ref}")	
		if(params["rc"]):
			shell("samtools faidx {input.query} {params.qrgns} | seqtk seq -r - > {output.query}")	
		else:	
			shell("samtools faidx {input.query} {params.qrgns} > {output.query}")	

rule RepeatMasker:
	input:
		fasta = temp("temp/{SM}_{SEQ}.fasta"),
	output:
		out = temp("temp/{SM}_{SEQ}.fasta.out"),
		cat = temp("temp/{SM}_{SEQ}.fasta.cat"),
		tbl = temp("temp/{SM}_{SEQ}.fasta.tbl"),
		msk = temp("temp/{SM}_{SEQ}.fasta.masked"),
	threads: 8
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (2 * threads),
        hrs = 72
	shell:
		"""
		RepeatMasker \
			-e ncbi \
			-species human \
			-dir $(dirname {input.fasta}) \
			-pa {threads} \
			{input.fasta}
		"""

rule DupMasker:
	input:
		fasta = temp("temp/{SM}_{SEQ}.fasta"),
		out = rules.RepeatMasker.output.out,
	output:
		dups = temp("temp/{SM}_{SEQ}.fasta.duplicons"),
	threads: 8
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (2 * threads),
        hrs = 72
	shell:
		"""
		DupMaskerParallel \
			-pa {threads} \
			-engine ncbi \
			{input.fasta}
		"""


rule DupMaskerColor:
	input:
		dups = rules.DupMasker.output.dups,
	output:
		dupcolor = temp("temp/{SM}_{SEQ}.fasta.duplicons.extra"),
	threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (1 * threads),
        hrs = 72
	shell:
		"""
		{SDIR}/scripts/DupMask_parserV6.pl -i {input.dups} -E -o {output.dupcolor}
		"""


rule get_cds:
	input:
		fasta = get_ref,
		bed = get_genes,
		ref = rules.get_rgns.output.ref,
	output:
		fasta = temp("temp/{SM}.genes.fasta"),
		bed = temp("temp/{SM}.ref.genes.bed"),
		bed12 = temp("temp/{SM}.ref.genes.12.bed"),
		tmp = temp("temp/tmp.{SM}.ref.genes.bed"),
	threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (1 * threads),
        hrs = 72
	params:
		bed = get_ref_bed,
	run:
		# read in regular symbol names 
		convert = {}
		for line in open(f"{SDIR}/data/gene_conversion.txt").readlines():
			t = line.strip().split()
			if(len(t) > 1):
				convert[  t[1].split(".")[0] ] = t[0]
	
		# make a bed file with all the corrdiantes in the correct space
		rtn = ""
		shell("> {output.fasta}")
		for bed in params["bed"]:
			#shell('bedtools intersect -a {input.bed} -b <(printf "{bed}") | bedtools bed12tobed6 -i /dev/stdin > {output.tmp}')
			shell('bedtools intersect -f 1.0 -a {input.bed} -b <(printf "{bed}") > {output.tmp}')
			
			# get all the cds seqeunces for map to the query 
			#shell("bedtools getfasta -name -split -fi {input.fasta} -bed {output.tmp} >> {output.fasta}")
			
			chrm, start, end = bed.split(); start = int(start);
			name = f"{chrm}:{start}-{end}"
			for line in open(output["tmp"]).readlines():
				t = line.strip().split()
				t[0] = name
				t[1] = int(t[1]) - start	
				t[2] = int(t[2]) - start	
				gene_id = t[3].split(".")[0]
				if(gene_id in convert): 
					t[3] = convert[gene_id]
					rtn += (11*"{}\t" + "{}\n").format(*t)
		open(output["bed12"], "w+").write(rtn)
		shell("bedtools bed12tobed6 -i {output.bed12} > {output.bed}")
		
		# get all the cds seqeunces to map to the query 
		shell("bedtools getfasta -s -name -split -fi {input.ref} -bed {output.bed12} > {output.fasta}")


rule query_genes:
	input:
		cds = rules.get_cds.output.fasta,
		query = rules.get_rgns.output.query,
	output:
		bam = temp("temp/{SM}.query.genes.bam",
		bed = temp("temp/{SM}.query.genes.bed",
		bed12 = temp("temp/{SM}.query.genes.12.bed",
	threads: 8
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (4 * threads),
        hrs = 72
	run:
	shell:
		"""
		minimap2 -p 0.98 --eqx -ax splice -C5 -O6,24 -B4 -t {threads} --cap-sw-mem=64g {input.query} {input.cds} | samtools view -b - | samtools sort - > {output.bam}
		bedtools bamtobed -bed12 -i {output.bam} > {output.bed12}
		bedtools bed12tobed6 -i {output.bed12} > {output.bed}
		"""

rule minimap2:
	input:
		ref = rules.get_rgns.output.ref,
		query = rules.get_rgns.output.query,
	output:
		paf = temp("temp/{SM}_{SCORE}_aln.paf"),
	params:
		score = get_score,
	threads: 16
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (4 * threads),
        hrs = 72
	shell:
		"""
		# YOU HAVE TO INCLUDE --cs FOR MINIMIRO TO WORK
		minimap2 -x asm20 -r 200000 -s {params.score} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}
		"""


rule minimiro:
	input:
		paf = rules.minimap2.output.paf,
		rmout = expand("temp/{{SM}}_{SEQ}.fasta.out", SEQ=SEQS),
		dmout = expand("temp/{{SM}}_{SEQ}.fasta.duplicons.extra", SEQ=SEQS),
		genes = rules.get_cds.output.bed12,
		query_genes = rules.query_genes.output.bed12,
	output:
		ps	= temp("temp/{SM}_{SCORE}_aln.ps"),
		pdf	= "minimiro/{SM}_{SCORE}_aln.pdf",
	threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (1 * threads),
        hrs = 72
	shell:"""
{SDIR}/minimiro.py --paf {input.paf} \
	--rm {input.rmout} \
	--dm {input.dmout} \
	--bed {input.genes} {input.query_genes} \
	--bestn 1000 \
	-o {output.ps} && \
	ps2pdf {output.ps} {output.pdf}
"""

rule coverage:
	input:
		bam = get_bam,
	output:
		png	= "minimiro/{SM}_coverage.png",
	threads: 1
    resources:
        mem = lambda wildcards, attempt, threads: attempt * (1 * threads),
        hrs = 72
	params:
		rgn = get_query_rgns, 
	shell:
		"""
		NucFreq/NucPlot.py {input.bam} {output.png} --regions {params.rgn} --height 4 --width 16
		"""

