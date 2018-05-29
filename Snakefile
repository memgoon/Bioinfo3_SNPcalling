import glob, os

### ENV
OrgName = "1000Genome"
suffix = ".fastq.gz"
prefix = "1.Fastq_dump/"
refSuffix = ".fa"

### FILES
RefFile = "0.data/Reference/hg38"+refSuffix


### VAR
SampleList = glob.glob(prefix+"*_1*"+suffix)
SampleList = [a_sample[len(prefix):-len(suffix)-2] for a_sample in SampleList]


### MAKE DIR
#os.system("mkdir 2.Trimmomatic")
#os.system("mkdir 3.Bowtie2")
#os.system("mkdir 4.Picard")
#os.system("mkdir 4.Picard/RG_bam")
#os.system("mkdir 4.Picard/RG_DU_bam")
#os.system("mkdir 4.Picard/RG_DU_FIX_bam")
#os.system("mkdir 5.GATK")
#os.system("mkdir 5.GATK/HaplotypeCaller")
#os.system("mkdir 5.GATK/RG_REaln_bam")
#os.system("mkdir 5.GATK/RG_REaln_REcal_bam")
#os.system("mkdir 5.GATK/SelectVariants")
#os.system("mkdir 6.VCF")

### ANNOUNCE
#print("Start run. Samples:")
#print(SampleList)

### RUN
rule all:
#	input: "5.GATK/VariantFiltration/{}.snp.vcf".format(OrgName), "5.GATK/VariantFiltration/{}.indel.vcf".format(OrgName)
	input: "6.VCF/final.vcf"

rule trimmomatic:
    input: prefix+"{GBSsample}"+suffix
    output: trimmed_fastq="2.Trimmomatic/{GBSsample}.fastq.gz", trimmed_log="2.Trimmomatic/{GBSsample}.log"
    params: adapter_type="TruSeq3-PE.fa"
    threads: 2
    shell: "java -jar /home/Program/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads {threads} -phred33 {input} {output.trimmed_fastq} ILLUMINACLIP:Program/Trimmomatic-0.36/adapters/{params.adapter_type}:2:30:10 TRAILING:20 MINLEN:80 2> {output.trimmed_log}"

rule makePicardDict:
	input: RefFile
	output: RefFile.replace(refSuffix, ".dict")
	shell: "java -Xmx8g -Djava.io.tmpdir=datatemp/ -jar /home/Program/Picard_2.8.1/picard.jar CreateSequenceDictionary REFERENCE={input} OUTPUT={output}"

rule makeBowtie2idx:
	input: RefFile
	params: RefFile.replace(refSuffix, "")
	output: RefFile.replace(refSuffix, ".rev.2.bt2")
	threads: 4
	shell: "bowtie2-build --threads {threads} {input} {params}"

rule makeFASTAidx:
	input: RefFile
	output: RefFile+".fai"
	shell: "samtools faidx {input}"

rule runBowtie2:
	input: fw_sample="2.Trimmomatic/{GBSsample}_1.fastq.gz", rv_sample="2.Trimmomatic/{GBSsample}_2.fastq.gz", reference=RefFile, ref_dict=RefFile.replace(refSuffix, ".dict"), ref_Bowtie2idx=RefFile.replace(refSuffix, ".rev.2.bt2"), ref_FASTAidx=RefFile+".fai"
	output: bam="3.Bowtie2/{GBSsample}.bam"
	params: ref_Bowtie2idx=RefFile.replace(refSuffix, "")
	threads: 4
	shell: "bowtie2 -p {threads} --no-unal -x {params.ref_Bowtie2idx} -1 {input.fw_sample} -2 {input.rv_sample}| samtools view -Sb - 1> {output.bam}"

rule getBowtie2stats:
    input: bam="3.Bowtie2/{GBSsample}.bam"
    output: stats="3.Bowtie2/{GBSsample}.stats"
    shell: "samtools stats {input.bam} > {output.stats}"

rule addReadGroup:
	input: bam="3.Bowtie2/{GBSsample}.bam", bwa_stats = ["3.Bowtie2/{}.stats".format(a_sample) for a_sample in SampleList]
	params: RGID="{GBSsample}", RGLB="{GBSsample}", RGSM="{GBSsample}"
	output: RG_bam="4.Picard/RG_bam/{GBSsample}.bam", RG_log="4.Picard/RG_bam/{GBSsample}.log"
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx8g -jar /home/Program/Picard_2.8.1/picard.jar AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.RG_bam} " \
	"SORT_ORDER=coordinate RGID={params.RGID} RGLB={params.RGLB} RGPL=illumina RGPU=non RGSM={params.RGSM} VALIDATION_STRINGENCY=LENIENT 2> {output.RG_log}"

rule removeDuplicates:
	input: RG_bam="4.Picard/RG_bam/{GBSsample}.bam"
	threads: 8 # To restrict too many jobs
	params: max_memory="48"
	output: RG_bam="4.Picard/RG_DU_bam/{GBSsample}.bam", metrics_filename="4.Picard/RG_DU_bam/{GBSsample}.metrics", RG_log="4.Picard/RG_DU_bam/{GBSsample}.log"
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx{params.max_memory}g -jar /home/Program/Picard_2.8.1/picard.jar MarkDuplicates INPUT={input.RG_bam} OUTPUT={output.RG_bam} " \
	"METRICS_FILE={output.metrics_filename} REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT 2> {output.RG_log}"

rule fixMate:
	input: RG_bam="4.Picard/RG_DU_bam/{GBSsample}.bam" 
	output: RG_FIX_bam="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam", RG_FIX_log="4.Picard/RG_DU_FIX_bam/{GBSsample}.log"
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx16g -jar /home/Program/Picard_2.8.1/picard.jar FixMateInformation INPUT={input.RG_bam} OUTPUT={output.RG_FIX_bam} " \
	"SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> {output.RG_FIX_log}"

rule indexBam_1:
	input: RG_bam="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam"
	output: RG_bai="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam.bai"
	shell: "samtools index {input.RG_bam}"

rule runRealign_TargetCreator:
	input: reference=RefFile, ref_dict=RefFile.replace(refSuffix, ".dict"),\
					 ref_idx=RefFile+".fai", RG_bam="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam",\
					 RG_bai="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam.bai"
	output: targetIntv = "5.GATK/RG_REaln_bam/{GBSsample}.intervals", REaln_log = "5.GATK/RG_REaln_bam/{GBSsample}.intervals.log"
	threads: 4 # to restrict
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx32g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T RealignerTargetCreator " \
	"-R {input.reference} -I {input.RG_bam} -o {output.targetIntv} 2> {output.REaln_log}"

rule runRealign_IndelRealigner:
	input: reference=RefFile, ref_idx=RefFile+".fai", RG_bam="4.Picard/RG_DU_FIX_bam/{GBSsample}.bam", targetIntv = "5.GATK/RG_REaln_bam/{GBSsample}.intervals"
	output: RG_REaln_bam="5.GATK/RG_REaln_bam/{GBSsample}.bam", REaln_log = "5.GATK/RG_REaln_bam/{GBSsample}.indelaligner.log"
	threads: 4 # to restrict
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx32g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T IndelRealigner " \
	"-R {input.reference} -I {input.RG_bam} -o {output.RG_REaln_bam} -targetIntervals {input.targetIntv} 2> {output.REaln_log}"

rule indexBam_2:
	input: realn_bam="5.GATK/RG_REaln_bam/{GBSsample}.bam"
	output: realn_bai="5.GATK/RG_REaln_bam/{GBSsample}.bam.bai"
	shell: "samtools index {input.realn_bam}"

rule baseQualityRecal:
	input: RG_REaln_bam="5.GATK/RG_REaln_bam/{GBSsample}.bam", RG_REaln_bai="5.GATK/RG_REaln_bam/{GBSsample}.bam.bai", reference=RefFile, vcf_Reference="0.data/dbSNP150/Homo_sapiens.vcf"
	output: recal_GRP="5.GATK/RG_REaln_REcal_bam/{GBSsample}.grp", recal_GRP_log="5.GATK/RG_REaln_REcal_bam/{GBSsample}.GRP_log"
	threads: 4 # to restrict
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx24g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T BaseRecalibrator " \
	"-R {input.reference} -I {input.RG_REaln_bam} -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate " \
	"-o {output.recal_GRP} -knownSites {input.vcf_Reference} 2> {output.recal_GRP_log}"

rule printReads:
	input: reference=RefFile, RG_REaln_bam="5.GATK/RG_REaln_bam/{GBSsample}.bam", recal_GRP="5.GATK/RG_REaln_REcal_bam/{GBSsample}.grp"
	output: recal_bam="5.GATK/RG_REaln_REcal_bam/{GBSsample}.bam", recal_BAM_log="5.GATK/RG_REaln_REcal_bam/{GBSsample}.bam.printReads_log"
	threads: 4 # to restrict
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx24g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T PrintReads " \
	"-R {input.reference} -I {input.RG_REaln_bam} -BQSR {input.recal_GRP} -o {output.recal_bam} 2> {output.recal_BAM_log}"

rule indexBam_3:
	input: realn_bam="5.GATK/RG_REaln_REcal_bam/{GBSsample}.bam"
	output: realn_bai="5.GATK/RG_REaln_REcal_bam/{GBSsample}.bam.bai"
	shell: "samtools index {input.realn_bam}"

### multisample vcf calling
rule haploTypeCaller_1:
	input: reference=RefFile, realn_bam=["5.GATK/RG_REaln_REcal_bam/{}.bam".format(a_sample) for a_sample in SampleList], \
	realn_bai=["5.GATK/RG_REaln_REcal_bam/{}.bam.bai".format(a_sample) for a_sample in SampleList]
	params: realn_bam=["-I 5.GATK/RG_REaln_REcal_bam/{}.bam".format(a_sample) for a_sample in SampleList]
	threads: 8
	output: raw_vcf="5.GATK/HaplotypeCaller/{}.vcf".format(OrgName), raw_vcf_log="5.GATK/HaplotypeCaller/{}.vcf.log".format(OrgName)
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx32g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R {input.reference} " \
	"{params.realn_bam} -o {output.raw_vcf} -nct {threads} 2> {output.raw_vcf_log}"

rule selectSNP:
	input: vcf="5.GATK/HaplotypeCaller/{}.vcf".format(OrgName), reference=RefFile
	output: snp="5.GATK/SelectVariants/{}.snp.vcf".format(OrgName), log="5.GATK/SelectVariants/{}.snp.log".format(OrgName)
	threads: 8
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx64g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T SelectVariants -R {input.reference} --variant {input.vcf} -o {output.snp} -selectType SNP -nt {threads} 2>> {output.log}"

rule selectINDEL:
	input: vcf="5.GATK/HaplotypeCaller/{}.vcf".format(OrgName), reference=RefFile
	output: indel="5.GATK/SelectVariants/{}.indel.vcf".format(OrgName), log="5.GATK/SelectVariants/{}.indel.log".format(OrgName)
	threads: 8
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx64g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -T SelectVariants -R {input.reference} --variant {input.vcf} -o {output.indel} -selectType INDEL -nt {threads} 2>> {output.log}"

rule INDELfilter:
	input: indel="5.GATK/SelectVariants/{}.indel.vcf".format(OrgName), reference = RefFile
	output: indel_filtered="5.GATK/VariantFiltration/{}.indel.vcf".format(OrgName), log="5.GATK/VariantFiltration/{}.indel.log".format(OrgName)
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx64g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -R {input.reference} -T VariantFiltration --variant {input.indel} -o {output.indel_filtered} " \
	"--filterExpression \"QUAL < 30\" --filterName \"QualFilter\" --filterExpression \"QD < 5.0\"  --filterName \"QD5\" --filterExpression \"FS > 200.0\" --filterName \"FS200\" 2> {output.log}"

rule SNPfilter:
	input: indel_filtered="5.GATK/VariantFiltration/{}.indel.vcf".format(OrgName), snp="5.GATK/SelectVariants/{}.snp.vcf".format(OrgName), reference = RefFile
	output: snp_filtered="5.GATK/VariantFiltration/{}.snp.vcf".format(OrgName), log="5.GATK/VariantFiltration/{}.snp.log".format(OrgName)
	shell: "java -Djava.io.tmpdir=datatemp/ -Xmx64g -jar /home/Program/GATK_3.7/GenomeAnalysisTK.jar -R {input.reference} -T VariantFiltration --variant {input.snp} -o {output.snp_filtered} --clusterSize 3 --clusterWindowSize 10 " \
	"--mask {input.indel_filtered} --maskName \"InDel\" --filterExpression \"QUAL < 30\" --filterName \"QualFilter\" --filterExpression \"FS > 200.0 \" --filterName \"FS200\" --filterExpression \"QD < 5.0\" --filterName \"QD5\" 2> {output.log}"

rule getPassSNP:
	input: snp_input="5.GATK/VariantFiltration/{}.snp.vcf".format(OrgName)
	output: snp_filtered="6.VCF/pass_snp_{}.vcf".format(OrgName)
	shell: "grep \"\#\|PASS\" {input} > {output}"

rule FilterVCFtools:
	input: "6.VCF/pass_snp_{}.vcf".format(OrgName)
	output: "6.VCF/final.recode.vcf"
	params: maf="0.05", miss="0.95", minMNDP="5"
	shell: "vcftools --vcf {input} --maf {params.maf} --max-missing {params.miss} --min-meanDP {params.minMNDP} --recode --out 6.VCF/final"

rule renameFinalVCF:
	input: snp_input="6.VCF/final.recode.vcf"
	output: snp_filtered="6.VCF/final.vcf"
	shell: "mv {input} {output}"
