configfile: 'config.yaml'
#	module load bedtools samtools sratoolkit/2.9.6 python/3.6.6 blast spades

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}


sra_by_name = {}
sra_by_path = {}

for s in sample_by_name.keys():
	if "SRA" in sample_by_name[s].keys():
		sra_by_name[s] = sample_by_name[s]["SRA"]
		sra_by_path[sample_by_name[s]["path"]] = sample_by_name[s]["SRA"]

sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]



def return_filename_by_sampname(sampname):
	filenames = []
	if sample_by_name[sampname]['paired']:
		filenames.append(sample_by_name[sampname]['readsfile1'])
		filenames.append(sample_by_name[sampname]['readsfile2'])
	else:
		filenames.append(sample_by_name[sampname]['readsfile'])
	return filenames


def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([pathprefix, fq]) for fq in filesin]
	return pathsout



# rule all:
# 	input: 
# 		pdf_out="results/1kGenHumanDenovoProjectOfDoom.pdf",
# 	params:
# 		runmem_gb=1,
# 		runtime="0:01:00",
# 		cores=1,
# 	run:
# 		shell(""" mkdir -p results/figures/; touch results/figures/null.png; for fig in results/figures/*png; do mv $fig $(echo $fig| rev | cut -f 2- -d . | rev ).$(date +%d_%b_%Y).png; done;  rm results/figures/null.*.png; """)
# 		shell(""" mkdir -p results/figures/supp/ ; touch results/figures/supp/null.png; for fig in results/figures/supp/*png; do mv $fig $(echo $fig| rev | cut -f 2- -d . | rev ).$(date +%d_%b_%Y).png; done; rm results/figures/supp/null.*.png; """)

# 		shell(""" mkdir -p results/tables/ ; touch results/tables/null.tmp ; for phial in $(ls -p results/tables/ | grep -v /); do pre=$(echo $phial | rev | cut -f 2- -d . | rev ); suff=$(echo $phial | rev | cut -f 1 -d . | rev ); mv results/tables/$phial results/tables/$pre.$(date +%d_%b_%Y).$suff; done ; rm results/tables/null.*.tmp; """)
# 		shell(""" mkdir -p results/tables/supp/ ; touch results/tables/supp/null.tmp ; for phial in $(ls -p results/tables/supp/ | grep -v /); do pre=$(echo $phial | rev | cut -f 2- -d . | rev ); suff=$(echo $phial | rev | cut -f 1 -d . | rev ); mv results/tables/supp/$phial results/tables/supp/$pre.$(date +%d_%b_%Y).$suff; done ; rm results/tables/supp/null.*.tmp; """)

# 		shell(""" mv results/1kGenHumanDenovoProjectOfDoom.pdf results/1kGenHumanDenovoProjectOfDoom.$(date +%d_%b_%Y).pdf """)
# 		shell(""" tar cf results.$(date +%d_%b_%Y).tar results/ """)


##########################	SUMMARIZE INPUT DATA 	###########################

rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "summaries/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
		clust_misc="",
	shell:
		"""
		mkdir -p summaries/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule summon_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("summaries/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "summaries/reference_genomes.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
		clust_misc="",
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"


##########################	GATHER SEQUENCED READS 	###########################

rule summon_reads_SRA_pe:
	output: 
		reads1='FASTQ/{path}/{prefix}_1.fastq',
		reads2='FASTQ/{path}/{prefix}_2.fastq',
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		clust_misc="",
	run:
		try:
			sra = sra_by_path[ "FASTQ/%s/" % tuple([wildcards.path]) ]
			shell(""" mkdir -p FASTQ/{wildcards.path}/ """)
#			shell("""
#				fasterq-dump --progress --split-3 --outdir FASTQ/{wildcards.path}/ --outfile {wildcards.prefix} {sra}
#			""")
			shell("""
				fastq-dump --defline-seq '@$sn[_$rn]/$ri' --outdir FASTQ/{wildcards.path}/ --split-3 {sra} 
			""")
		except KeyError:
			raise KeyError("Sample is listed as empirical but no reads found and no SRA to download!" )

rule summon_reads_SRA_se:
	output: 
		reads='FASTQ/{path}/{prefix,SRR\d+}.fastq',
	params:
		runmem_gb=16,
		runtime="3:00:00",
		cores=1,
		clust_misc="",
#	wildcard_constraints:
#        prefix="SRR\d+",
	run:
		try:
			sra = sra_by_path[ "FASTQ/%s/" % tuple([wildcards.path]) ]
			shell(""" mkdir -p FASTQ/{wildcards.path}/ """)
#			shell("""
#				fasterq-dump --progress --split-3 --outdir FASTQ/{wildcards.path}/ --outfile {wildcards.prefix} {sra}
#			""")
			shell("""
				fastq-dump --defline-seq '@$sn[_$rn]/$ri' --outdir FASTQ/{wildcards.path}/ --split-3 {sra} 
			""")
		except KeyError:
			raise KeyError("Sample is listed as empirical but no reads found and no SRA to download!" )


# rule fastp_clean_sample_pe:
# 	input:
# 		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
# 	output:
# 		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
# 		jason = "{pathprefix}/{samplename}.json"
# 	params:
# 		runmem_gb=8,
# 		runtime="6:00:00",
# 		cores=1,
# 		clust_misc="",
# 		#--trim_front1 and -t, --trim_tail1
# 		#--trim_front2 and -T, --trim_tail2. 
# 		common_params = "--json {pathprefix}/{samplename}.json",# --html meta/FASTP/{samplename}.html", 
# 		pe_params = "--detect_adapter_for_pe --correction",
# 	message:
# 		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
# 	shell:
# 		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"



rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "{pathprefix}/{samplename}.False.json"
	params:
		runmem_gb=16,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.False.json",# --html meta/FASTP/{samplename}.html", 
		se_params = "",
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "{pathprefix}/{samplename}.True.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.True.json",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"









# rule FASTP_summarizer:
# 	input: 
# 		jason = lambda wildcards: expand("{path}{samp}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, )
# 	output:
# 		jason_pruned = "summaries/FASTP/{samplename}.json.pruned"
# 	params:
# 		runmem_gb=1,
# 		runtime="5:00",
# 		cores=1,
# 		clust_misc="",
# 	message:
# 		"Summarizing reads for sample ({wildcards.samplename}) .... "	
# 	shell:
# 		"""
# 		cp {input.jason} summaries/FASTP/{wildcards.samplename}.json
# 		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
# 		"""

rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}{samp}.{pairt}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
	output:
		jason_pruned = "summaries/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		mkdir -p summaries/FASTP/
		cp {input.jason} summaries/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""


rule summon_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("summaries/FASTP/{samplename}.json.pruned", samplename=sampname_by_group['all'])
	output:
		summary = "summaries/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
		clust_misc="",
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"

def read_numbs(twoo):
	if twoo:
		return([1,2])
	else:
		return([0])

rule spades_assembler:
	input:
		reads_in = lambda wildcards: expand("{path}{sample}.clean.R{pair}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, pair = read_numbs(sample_by_name[wildcards.sample]['paired']) ),
	output:
		scaff = "genome_assemblies/spades/{sample}/{sample}.spades.scaffolds.fa"
	params:
		runmem_gb=512,
		runtime="196:00:00",
		cores=1,
		cpus_per_task=12,
		clust_misc="--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" mkdir -p genome_assemblies/spades/{wildcards.sample}/""")

		if sample_by_name[wildcards.sample]['paired']:
			rudd = " -1 %s -2 %s " %tuple([input.reads_in[0], input.reads_in[1]])
		else:
			rudd = " -s %s " %tuple([input.reads_in[0]])

		shell("""spades.py -o genome_assemblies/spades/{wildcards.sample}/ {rudd} -t {params.cpus_per_task} -m {params.runmem_gb}""")
#sbatch --mem=196G --partition=bigmem --qos bigmem_access --out=results/xpodTest.illumina_spades/spades.out --ntasks=1 --cpus-per-task=8 -N 1 --time=96:00:00 --wrap="spades.py -o results/xpodTest.illumina_spades -1 /proj/cdjones_lab/csoeder/xpod/fastq/xLae/SRR3210978_1.fastq -2 /proj/cdjones_lab/csoeder/xpod/fastq/xLae/SRR3210978_2.fastq  -t 8 -m 196"
		shell(""" mv genome_assemblies/spades/{wildcards.sample}/scaffolds.fasta {output.scaff} """)


# use quast to analyze 



rule make_blastDB_from_genome:
	input:
		genome_in = "genome_assemblies/spades/{sample}/{sample}.spades.scaffolds.fa"
	output:
		blastdb_flag = "utils/flags/blastdb_made/genome_assemblies/{sample}/{sample}.spades.scaffolds.blastdbmade.flag",
	params:
		runmem_gb=8,
		runtime="45:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" mkdir -p utils/flags/blastdb_made/genome_assemblies/{wildcards.sample}/ """)
		shell(""" makeblastdb -in {input.genome_in} -dbtype nucl """)
		shell(""" touch {output.blastdb_flag}""")



rule scan_assemblies_for_xpod:
	input:
		genome_in = "genome_assemblies/spades/{sample}/{sample}.spades.scaffolds.fa",
		blastdb_flag = "utils/flags/blastdb_made/genome_assemblies/{sample}/{sample}.spades.scaffolds.blastdbmade.flag"

	output:
		pept_out ="genome_assemblies/spades/{sample}/{sample}.peptBLAST.out",
		nucl_out = "genome_assemblies/spades/{sample}/{sample}.nuclBLAST.out"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell("""tblastn -query utils/xpod.pept.fa -db {input.genome_in} -outfmt 6 > {output.pept_out} """)
		shell("""blastn -query utils/xpod.nucl.fa -db {input.genome_in} -outfmt 6 > {output.nucl_out} """)


rule scan_all_assemblies:
	input:
		peppys = expand("genome_assemblies/spades/{sample}/{sample}.peptBLAST.out", sample=sampname_by_group['all']),
		nukular = expand("genome_assemblies/spades/{sample}/{sample}.nuclBLAST.out", sample=sampname_by_group['all']),
	output:
		flagout = "utils/flags/genomes_scanned_for_xpod.flag"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" mkdir -p utils/flags/ """)
		shell(""" touch {output.flagout} """)





# rule write_report:
# 	input:
# 		reference_genome_summary = ["summaries/reference_genomes.summary"],
# 		reference_annotation_summary = ["summaries/reference_annotations.summary"],
# 		gene_lists_summary = ["summaries/geneLists.stats"],
# 		sequenced_reads_summary=["summaries/sequenced_reads.dat"],
# 		aligned_reads_summary = expand("summaries/alignments.vs_dm6main.{aligner}.summary", aligner=["mapspliceRaw","mapspliceMulti","mapspliceUniq","mapspliceRando"]),
# 		counting_summary_genes = expand("summaries/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.stat", group ="all", ref_genome = "dm6main", annot = ["dm6_genes","fru_exons"], flag=["MpBC","MpBCO"], aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"]),#"mapspliceRaw",
# 		counting_summary_exons = expand("summaries/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.stat", group ="all", ref_genome = "dm6main", annot = ["fru_junct", "fru_intron"], flag=["MpBC","MpBCO"], aligner=["mapspliceMulti_SpliceOnly","mapspliceUniq_SpliceOnly","mapspliceRando_SpliceOnly"]),#"mapspliceRaw",
# 		expFlg = "utils/all.expression.flag",
# 		diff_exprs = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["grpWtVs47b","grpWtVs67d","grpWtVsFru","grpWtVsMut","hausWtVsMut","wildTypeHousing"], annot = ["dm6_genes","fru_exons"], flag= ["MpBC", "MpBCO"],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["de"]),#,"itemized.de"]),
# 		diff_item = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["hausWtVsMut", "hausWtVsMut_noFru"], annot = ["dm6_genes","fru_exons"], flag= ["MpBC", "MpBCO"],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["itemized.de"]),#,""]),
# 		diff_item2 = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["2_days_difference", "47b_on_88a", "cantonAmos", "cantonWt", "cantonAmosWt"], annot = ["dm6_genes","fru_exons"], flag= ["MpBC", ],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["itemized.de"]),#,""]),
# 		diff_fru = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["hausWtVsMut", "hausWtVsMut_noFru"], annot = ["fru_junct", "fru_intron", "fru_exons"], flag= [ "MpBCO"],aligner=["mapspliceMulti_SpliceOnly","mapspliceUniq_SpliceOnly","mapspliceRando_SpliceOnly"], suff = ["itemized.de"]),#,""]),
# 		gene_onts = expand("gene_ont/{contrast}/{contrast}.vs_dm6main.dm6_genes.mapsplice{aligner}.MpBC.go", aligner = ["Multi","Rando","Uniq"], contrast = ["hausWtVsMut","47b_on_88a","2_days_difference", "cantonAmos", "cantonWt", "cantonAmosWt"]),

# 	output:
# 		pdf_out="results/VolkanLab_BehaviorGenetics.pdf",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=2,
# 	message:
# 		"writing up the results.... "
# 	run:
# 		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
# 		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
# 		shell("""mkdir -p results/figures/supp/ results/tables/supp/""")
# 		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/RNAseq_results.Rmd',output_file='{pwd}{output.pdf_out}')"  """)
# #		shell(""" tar cf results.tar results/ """)


