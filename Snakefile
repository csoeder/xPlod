configfile: 'config.yaml'
#	module load bedtools samtools sratoolkit/2.9.6 python/3.6.6 blast spades quast r blat fastx_toolkit

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
ref_genome_by_species = { g['species'] : g for g in config['reference_genomes']}


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
		fa_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['path'],

	output:
		report_out = "summaries/reference_genomes/{ref_gen}.fai.report",
		quast_out = "summaries/reference_genomes/quast/{ref_gen}/report.tsv",
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
		clust_misc="",
		quast_params  = " --eukaryote ",
	shell:
		"""
		mkdir -p summaries/reference_genomes/quast/{wildcards.ref_gen}/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		quast.py  {params.quast_params} --output-dir summaries/reference_genomes/quast/{wildcards.ref_gen}/ {input.fa_in} 
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


##########################	SEARCH IN VAIN FOR THE SEQUENCE 	###########################
##########################	IN THE REFERENCE GENOMES	###########################

rule sourceESTblatterUp:
	input:
		ests = "utils/source_ESTs.fa",
	output:
		raw_blat = "scans/background/source_ESTs.blat.vs_{ref_gen}.mapt.psl",
		raw_covBed = "scans/background/source_ESTs.blat.vs_{ref_gen}.mapt.queryCoverage.bed",

		filtered_blat = "scans/background/source_ESTs.blat.vs_{ref_gen}.mapt.filt.psl",
		filtered_covbed = "scans/background/source_ESTs.blat.vs_{ref_gen}.mapt.filt.queryCoverage.bed",
	params:
		runmem_gb=16,
		runtime="24:00:00",
		cores=8,
		match_thresh = 50,#75 bases must match
		match_frac = 0.50,#75% of query must align
	message:
		"blatting source_ESTs against reference_genome {wildcards.ref_gen} .... "
	run:
		shell("""mkdir -p scans/background/""")
		ref_gen = ref_genome_by_name[wildcards.ref_gen]["path"]
		shell(""" blat -noHead {ref_gen} {input.ests} {output.raw_blat}; """)
		shell(""" cat {output.raw_blat}| awk '{{if(($1)>{params.match_thresh})print;}}' | awk '{{if(($13-$12)/$11>{params.match_frac})print;}}' |  > {output.filtered_blat} """)
#		shell(""" ~/modules/UCSC_utils/pslToBed {output.filtered_transcriptome_blat} {output.transcriptome_bed} """)
#		shell(""" cat {input.transcriptome_in} | fasta_formatter -t | grep -wFf <( cat {output.raw_transcriptome_blat} | cut -f 10 | sort | uniq ) | awk '{{print">"$0}}' | tr "\t" "\n" > {output.mappable_seqs}; """)
#		shell(""" cat {input.transcriptome_in} | fasta_formatter -t | grep -wFf <( cat {output.filtered_transcriptome_blat} | cut -f 10 | sort | uniq ) | awk '{{print">"$0}}' | tr "\t" "\n" > {output.filtered_mappable_seqs}; """)

		shell("""cat {output.raw_blat} | awk '{{print$10,$12,$13,$14,$1,$9}}' | tr " " "\t" > {output.raw_covBed}""")
		shell("""cat {output.filtered_blat} | awk '{{print$10,$12,$13,$14,$1,$9}}' | tr " " "\t" > {output.filtered_covbed}""")






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



rule assembly_quasticator:
	input: 
		assemble_in = "genome_assemblies/spades/{sample}/{sample}.spades.scaffolds.fa"
	output:
		quast_out = "genome_assemblies/spades/{sample}/quast/report.tsv"
	params:
		runmem_gb=64,
		runtime="32:00:00",
		cores=1,
		quast_params  = " --eukaryote  --fragmented --conserved-genes-finding --no-html --no-plots ",
	message:
		"Summarizing reads for sample ({wildcards.sample}) .... "	
	run:
		sample_stuff = sample_by_name[wildcards.sample]
		if sample_stuff["paired"]:
			read_str = " -1 %s%s -2 %s%s " % tuple([sample_stuff["path"],sample_stuff["readsfile1"],sample_stuff["path"],sample_stuff["readsfile2"]])
		elif sample_stuff["platform"] == "pacbio":
			read_str = " --pacbio %s%s " % tuple([sample_stuff["path"],sample_stuff["readsfile"]])
		elif sample_stuff["platform"] == "nanopore":
			read_str = " --nanopore %s%s " % tuple([sample_stuff["path"],sample_stuff["readsfile"]])
		else: 
			read_str = " --single %s%s " % tuple([sample_stuff["path"],sample_stuff["readsfile"]])
		ref_gen = ref_genome_by_species[sample_stuff["species"][0]]["path"]
		shell(""" quast.py -r {ref_gen} {params.quast_params} {read_str} --output-dir genome_assemblies/spades/{wildcards.sample}/quast/ {input.assemble_in}  """)


rule quasticate_everything:
	input: 
		quasts = lambda wildcards: expand( "genome_assemblies/spades/{sample}/quast/report.tsv", sample=sampname_by_group['all'] )
	output:
		quast_out = "summaries/genome_assembly.summary"
	params:
		runmem_gb=8,
		runtime="5:00",
		cores=1,
	message:
		"quasitcating....... .... "	
	run:
		for samp in sampname_by_group['all']:
			shell(""" cat genome_assemblies/spades/{sample}/quast/report.tsv | awk '{{print"{sample}\t"$0}}' > {output.quast_out} """)
#	add -r for genome reference??



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
		pept_out = "scans/assemblies/spades/{sample}.peptBLAST.out",
		nucl_out = "scans/assemblies/spades/{sample}.nuclBLAST.out",
		est_out = "scans/assemblies/spades/{sample}.estBLAST.out",

		nucl_qCov = "scans/assemblies/spades/{sample}.nuclBLAST.queryCoverage.bed",
		est_qCov = "scans/assemblies/spades/{sample}.estBLAST.queryCoverage.bed",

	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc="",#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
		evalue = 2,

	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" mkdir -p scans/assemblies/spades/ """)
		shell("""tblastn -query utils/xpod.pept.fa -db {input.genome_in} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -evalue {params.evalue} > {output.pept_out} """)
		
		shell("""blastn -query utils/xpod.nucl.fa -db {input.genome_in} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -evalue {params.evalue} > {output.nucl_out} """)
		shell("""cat {output.nucl_out} | awk '{{print$1,$7,$8,$2,$12,$13}}' | tr " " "\t" | sed -e 's/plus/+/g' | sed -e 's/minus/-/g' > {output.nucl_qCov}""")

		shell("""blastn -query utils/source_ESTs.fa -db {input.genome_in} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -evalue {params.evalue} > {output.est_out} """)
		shell("""cat {output.est_out} | awk '{{print$1,$7,$8,$2,$12,$13}}' | tr " " "\t" | sed -e 's/plus/+/g' | sed -e 's/minus/-/g' > {output.est_qCov}""")


# cat scans/assemblies/spades/xLae1.nuclBLAST.out | awk '{print$1,$7,$8,$2,$12,$13,"xLae1"}' | tr " " "\t" | sed -e 's/plus/+/g' | sed -e 's/minus/-/g' 


rule scan_all_assemblies:
	input:
		peppys = expand("scans/assemblies/spades/{sample}.peptBLAST.out", sample=sampname_by_group['all']),
		nukular = expand("scans/assemblies/spades/{sample}.nuclBLAST.out", sample=sampname_by_group['all']),
		ests = expand("scans/assemblies/spades/{sample}.estBLAST.out", sample=sampname_by_group['all']),

	output:
		statout = "summaries/scans/xpod_vs_assemblies.stat"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" mkdir -p summaries/scans/ ; rm -rf {output.statout} ; """)

		for samp in sampname_by_group['all']:
			species = sample_by_name[samp]["species"][0]
			shell(""" cat scans/assemblies/spades/{samp}.peptBLAST.out | wc -l | awk '{{print"{species}\t{samp}\traw\tpeptide\t"$0}}' >> {output.statout}  """)
			shell(""" cat scans/assemblies/spades/{samp}.nuclBLAST.out | wc -l | awk '{{print"{species}\t{samp}\traw\tnucleotide\t"$0}}' >> {output.statout}  """)
			shell(""" cat scans/assemblies/spades/{samp}.estBLAST.out | wc -l | awk '{{print"{species}\t{samp}\traw\tEST\t"$0}}' >> {output.statout}  """)




rule FOOG_extractor:#sort BLAST by query start on theory that exon bounaries will come out consistently
	input:
		nukular = "scans/assemblies/spades/{sample}.nuclBLAST.out",

	output:
		foogBed = "scans/assemblies/spades/{sample}.nuclBLAST.FOOG.bed",
		foogFa = "scans/assemblies/spades/{sample}.nuclBLAST.FOOG.fa",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		#shell(""" mkdir -p summaries/scans/ ; rm -rf {output.statout} ; """)
		shell("""cat {input.nukular}  | sort -k 7 -g | cut -f 2,9,10 | awk '{{print $1,$2<$3?$2:$3,$3<$2?$2:$3,"FOOG_"NR,"{wildcards.sample}"}}' | tr " " "\t"  > {output.foogBed} """)
		shell("""bedtools getfasta -fi genome_assemblies/spades/{wildcards.sample}/{wildcards.sample}.spades.scaffolds.fa -fo - -bed  {output.foogBed} -name | sed -e "s/>F/>{wildcards.sample}:F/g" > {output.foogFa} """)



rule FOOGazi_repeater:#sort BLAST by query start on theory that exon bounaries will come out consistently
	input:
		foogies = expand("scans/assemblies/spades/{sample}.nuclBLAST.FOOG.fa", sample=sampname_by_group['all']),
		beds = expand("scans/assemblies/spades/{sample}.nuclBLAST.FOOG.bed", sample=sampname_by_group['all']),

	output:
		flag = "utils/foog.flag",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
		clust_misc=""#--cpus-per-task=%s --partition=bigmem --qos bigmem_access" %tuple([8])#[params.cpus_per_task]),
	message:
		"Collecting read summaries for all samples ...."
	run:
		shell(""" cat {input.beds} | cut -f 4 | sort | uniq >  scans/assemblies/spades/FOOGs.list """)
		shell(""" mkdir -p scans/assemblies/spades/FOOGs/ """)
		shell("""cat scans/assemblies/spades/FOOGs.list | while read -r foog; do cat {input.foogies} | fasta_formatter -t | tr ":" "\t" | grep -w "$foog" | awk '{{print">"$1;print$3}}' > scans/assemblies/spades/FOOGs/"$foog".fa ; done""")
		shell( """ touch {output.flag}""")


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


