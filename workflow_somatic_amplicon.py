#!/usr/bin/env python

# Standard packages
import os
import sys
import getpass
import argparse
import configuration
import subprocess as sub

from toil.job import Job
from toil.job import JobException
from cassandra.auth import PlainTextAuthProvider


def run_and_log_command(command, logfile):
    """This function uses the python subprocess method to run the specified command and writes all error to the
    specified logfile

    :param command: The command-line command to execute.
    :type name: str.
    :param logfile: The logfile to output error messages to.
    :type logfile: str.
    :returns:  Nothing
    :raises: RuntimeError

    """

    with open(logfile, "wb") as err:
        sys.stdout.write("Executing {} and writing to logfile {}\n".format(command, logfile))
        err.write("Command: {}\n".format(command))
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            raise RuntimeError("An error occurred when executing the commandline: {}. "
                               "Please check the logfile {} for details\n".format(command, logfile))


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")


def run_fastqc(job, config, samples):
    """Run FastQC on provided FastQ files
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: Samples dictionary
    :type samples: str.
    """

    job.fileStore.logToMaster("Running FastQC for all samples\n")
    logfile = "fastqc.log"

    fastq_files_list = list()
    for sample in samples:
        fastq_files_list.append(samples[sample]['fastq1'])
        fastq_files_list.append(samples[sample]['fastq2'])

    fastq_files_string = " ".join(fastq_files_list)
    command = ["{}".format(config['fastqc']['bin']),
               "{}".format(fastq_files_string),
               "--extract"]

    job.fileStore.logToMaster("FastQC Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)


def run_bwa_mem(job, config, name, samples):
    """Run GATK's DiagnoseTargets against the supplied region

    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param fastq1: Input FastQ File.
    :type fastq1: str.
    :param fastq2: Input FastQ File.
    :type fastq2: str.
    :returns:  str -- Aligned and sorted BAM file name.

    """

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(name))

    output_bam = "{}.bwa.sorted.bam".format(name)
    temp = "{}.bwa.sort.temp".format(name)
    logfile = "{}.bwa-align.log".format(name)

    bwa_cmd = ["{}".format(config['bwa']['bin']),
               "mem",
               "-t",
               "{}".format(config['bwa']['num_cores']),
               "-M",
               "-v",
               "2",
               "{}".format(config['reference']),
               "{}".format(samples[name]['fastq1']),
               "{}".format(samples[name]['fastq2'])]

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(config['bwa']['num_cores']),
                "-O",
                "bam",
                "-o",
                "{}".format(output_bam),
                "-T",
                "{}".format(temp),
                "-"]

    command = "{} | {} | {}".format(" ".join(bwa_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("BWA Command: {}\n".format(command))
    run_and_log_command(command, logfile)

    return output_bam


def add_or_replace_readgroups(job, config, name, input_bam):
    """Run Picard's AddOrReplaceReadGroups on the specified BAM
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.
    """

    job.fileStore.logToMaster("Running AddOrReplaceReadGroups in sample: {}".format(name))

    output_bam = "{}.rg.sorted.bam".format(name)
    logfile = "{}.addreadgroups.log".format(name)
    index_log = "{}.buildindex.log".format(name)

    command = ["{}".format(config['picard-add']['bin']),
               "AddOrReplaceReadGroups",
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam),
               "RGID={}".format(name),
               "RGSM={}".format(name),
               "RGLB={}".format(name),
               "RGPL=illumina",
               "RGPU=miseq"]

    command2 = ["{}".format(config['picard-add']['bin']),
                "BuildBamIndex",
                "INPUT={}".format(output_bam)]

    job.fileStore.logToMaster("GATK AddOrReplaceReadGroupsCommand Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command2))
    run_and_log_command(" ".join(command2), index_log)

    return output_bam


def realign_target_creator(job, config, name, input_bam):
    """Run GATK TargetCreator on the specified BAM to identify targets for realignment
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The file name of the targets file.
    """

    targets = "{}.targets.intervals".format(name)
    targets_log = "{}.targetcreation.log".format(name)

    command = ["{}".format(config['gatk-realign']['bin']),
               "-T",
               "RealignerTargetCreator",
               "-R",
               "{}".format(config['reference']),
               "-I",
               "{}".format(input_bam),
               "-o",
               "{}".format(targets),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-nt",
               "{}".format(config['gatk-realign']['num_cores'])
               ]

    job.fileStore.logToMaster("GATK RealignerTargetCreator Command: {}\n".format(command))
    run_and_log_command(" ".join(command), targets_log)

    return targets


def realign_indels(job, config, name, input_bam, targets):
    """Run GATK Indel Realignment on the specified BAM
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :param targets: The file name of targets to realign.
    :type targets: str.
    :returns:  str -- The output bam file name.
    """

    output_bam = "{}.realigned.sorted.bam".format(name)
    realign_log = "{}.realignindels.log".format(name)

    command = ["{}".format(config['gatk-realign']['bin']),
               "-T",
               "IndelRealigner",
               "-R",
               "{}".format(config['reference']),
               "-I",
               "{}".format(input_bam),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-targetIntervals",
               "{}".format(targets),
               "--read_filter",
               "NotPrimaryAlignment",
               "-o",
               "{}".format(output_bam)]

    job.fileStore.logToMaster("GATK IndelRealigner Command: {}\n".format(command))
    run_and_log_command(" ".join(command), realign_log)

    return output_bam


def recalibrator(job, config, name, input_bam):
    """Run GATK Recalibrator on the specified BAM

    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.

    """

    output_bam = "{}.recalibrated.sorted.bam".format(name)
    recal_config = "{}.recal".format(name)
    recal_log = "{}.recalibrate.log".format(name)
    print_log = "{}.printrecalibrated.log".format(name)
    cp_log = "{}.copy.log".format(name)

    # Calculate covariates
    recal_commands = ["java -Xmx8g -jar /mnt/shared-data/anaconda2/envs/ddb/bin/GenomeAnalysisTK.jar",
                      #"{}".format(config['gatk-recal']['bin']),
                      "-T",
                      "BaseRecalibrator",
                      "-R",
                      "{}".format(config['reference']),
                      "-I",
                      "{}".format(input_bam),
                      "-o",
                      "{}".format(recal_config),
                      "--knownSites",
                      "{}".format(config['dbsnp']),
                      "-nct",
                      "{}".format(config['gatk-recal']['num_cores'])]

    # Print recalibrated BAM
    print_reads_command = ["{}".format(config['gatk-recal']['bin']),
                           "-T",
                           "PrintReads",
                           "-R",
                           "{}".format(config['reference']),
                           "-I",
                           "{}".format(input_bam),
                           "-o",
                           "{}".format(output_bam),
                           "-BQSR",
                           "{}".format(recal_config),
                           "-nct",
                           "{}".format(config['gatk-recal']['num_cores'])]

    # Copy index to alternative name
    cp_command = ["cp",
                  "{}.recalibrated.sorted.bai".format(name),
                  "{}.recalibrated.sorted.bam.bai".format(name)]

    job.fileStore.logToMaster("GATK BaseRecalibrator Command: {}\n".format(recal_commands))
    run_and_log_command(" ".join(recal_commands), recal_log)

    job.fileStore.logToMaster("GATK PrintReads Command: {}\n".format(print_reads_command))
    run_and_log_command(" ".join(print_reads_command), print_log)

    job.fileStore.logToMaster("GATK Copy Command: {}\n".format(cp_command))
    run_and_log_command(" ".join(cp_command), cp_log)

    return output_bam


def sambamba_region_coverage(job, config, name, samples, input_bam):
    """Run SamBambam to calculate the coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample/library name.
    :type name: str.
    :param input_bam: The input_bam file name to process.
    :type samples: dict
    :param samples: The samples configuration dictionary
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.sambamba_coverage.bed".format(name)
    logfile = "{}.sambamba_coverage.log".format(name)

    command = ["{}".format(config['sambamba']['bin']),
               "depth region",
               "-L",
               "{}".format(samples[name]['regions']),
               "-t",
               "{}".format(config['sambamba']['num_cores']),
               "-T",
               "{}".format(config['coverage_threshold']),
               "-T",
               "{}".format(config['coverage_threshold2']),
               "{}".format(input_bam),
               ">",
               "{}".format(output)]

    job.fileStore.logToMaster("SamBamba Coverage Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)

    return output


def freebayes_single(job, config, name, input_bam):
    """Run FreeBayes without a matched normal sample
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    freebayes_vcf = "{}.freebayes.vcf".format(name)
    logfile = "{}.freebayes.log".format(name)

    command = ["{}".format(config['freebayes']['bin']),
               "--fasta-reference",
               "{}".format(config['reference']),
               "--min-alternate-fraction",
               "{}".format(config['min_alt_af']),
               "--pooled-discrete",
               "--pooled-continuous",
               "--genotype-qualities",
               "--report-genotype-likelihood-max",
               "--allele-balance-priors-off",
               "--use-duplicate-reads",
               "--min-repeat-entropy 1",
               "-v",
               "{}".format(freebayes_vcf),
               "{}".format(input_bam)]

    job.fileStore.logToMaster("FreeBayes Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)

    return freebayes_vcf


def mutect_single(job, config, name, samples, input_bam):
    """Run MuTect on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param samples: samples configuration dictionary
    :type samples: dict
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    mutect_vcf = "{}.mutect.vcf".format(name)
    temp_mutect = "{}.tempmutect.vcf".format(name)

    output_stats = "{}.mutectstats.txt".format(name)
    sample_coverage = "{}.mutectcoverage.wig.txt".format(name)

    mutect_logfile = "{}.mutect.log".format(name)
    subset_log = "{}.mutect_subset.log".format(name)

    mutect_command = ["{}".format(config['mutect']['bin']),
                      "-T",
                      "MuTect",
                      "-R",
                      "{}".format(config['reference']),
                      "--dbsnp",
                      "{}".format(config['dbsnp']),
                      "--cosmic",
                      "{}".format(config['cosmic']),
                      "--enable_extended_output",
                      "-I:tumor",
                      "{}".format(input_bam),
                      "--coverage_file",
                      "{}".format(sample_coverage),
                      "-L",
                      "{}".format(samples[name]['regions']),
                      "-isr",
                      "INTERSECTION",
                      "-im",
                      "ALL",
                      "-dt",
                      "NONE",
                      "-o",
                      "{}".format(output_stats),
                      "-vcf",
                      "{}".format(temp_mutect)]

    subset_command = ["cat",
                      "{}".format(temp_mutect),
                      "|",
                      "{}".format(config['vcftools_subset']['bin']),
                      "-e",
                      "-c",
                      "{}".format(name),
                      ">",
                      "{}".format(mutect_vcf)]

    job.fileStore.logToMaster("MuTect Command: {}\n".format(mutect_command))
    run_and_log_command(" ".join(mutect_command), mutect_logfile)

    job.fileStore.logToMaster("Subset Command: {}\n".format(subset_command))
    run_and_log_command(" ".join(subset_command), subset_log)

    return mutect_vcf


def platypus_single(job, config, name, samples, input_bam):
    """Run Platypus on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    platypus_vcf = "{}.platypus.vcf".format(name)
    platypus_log = "{}.platypus.log".format(name)
    internal_log = "{}.platypus_internal.log".format(name)

    platypus_command = ["{}".format(config['platypus']['bin']),
                        "callVariants",
                        "--refFile={}".format(config['reference']),
                        "--regions={}".format(samples[name]['regions']),
                        "--assemble=1",
                        "--assembleBadReads=1",
                        "--assembleBrokenPairs=1",
                        "--filterDuplicates=0",
                        "--minVarFreq={}".format(config['min_alt_af']),
                        "--nCPU={}".format(config['platypus']['num_cores']),
                        "--logFileName={}".format(internal_log),
                        "--bamFiles={}".format(input_bam),
                        "--output={}".format(platypus_vcf)]

    job.fileStore.logToMaster("Platypus Command: {}\n".format(platypus_command))
    run_and_log_command(" ".join(platypus_command), platypus_log)

    return platypus_vcf


def scalpel_single(job, config, name, samples, input_bam):
    """Run Scalpel on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    cwd = os.getcwd()
    output_dir = os.path.join(cwd, "{}-scalpel-output".format(name))
    scalpel_vcf = os.path.join(output_dir, "variants.indel.vcf")
    fixed_vcf = "{}.scalpel.vcf".format(name)
    logfile = "{}.scalpel.log".format(name)
    logfile2 = "{}.scalpel_fix.log".format(name)

    scalpel_command = ["{}".format(config['scalpel']['bin']),
                       "--single",
                       "--intarget",
                       # "--covthr",
                       # "3",
                       # "--lowcov",
                       # "1",
                       "--ref",
                       "{}".format(config['reference']),
                       "--bed",
                       "{}".format(samples[name]['regions']),
                       "--format",
                       "vcf",
                       "--numprocs",
                       "{}".format(config['scalpel']['num_cores']),
                       "--bam",
                       "{}".format(input_bam),
                       "--dir",
                       "{}".format(output_dir)]

    fix_sample_name_command = ["cat",
                               "{}".format(scalpel_vcf),
                               "|",
                               "sed",
                               "'s/sample/{}/g'".format(name),
                               ">",
                               "{}".format(fixed_vcf)]

    job.fileStore.logToMaster("Scalpel Command: {}\n".format(scalpel_command))
    run_and_log_command(" ".join(scalpel_command), logfile)

    job.fileStore.logToMaster("Scalpel Fix Command: {}\n".format(fix_sample_name_command))
    run_and_log_command(" ".join(fix_sample_name_command), logfile2)

    file_path = os.path.join(cwd, fixed_vcf)
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return scalpel_vcf
    else:
        job.fileStore.logToMaster("Scalpel ran into a problem and no output was generated for file {}. Check logfile"
                                  "{} for details\n".format(scalpel_vcf, logfile))
        return JobException("Scalpel ran into a problem and no output was generated for file {}. Check logfile"
                            "{} for details\n".format(scalpel_vcf, logfile))


def vardict_single(job, config, name, samples, input_bam):
    """Run VarDict on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    vardict_vcf = "{}.vardict.vcf".format(name)
    logfile = "{}.vardict.log".format(name)

    vardict = ["{}".format(config['vardict']['bin']),
               "-G",
               "{}".format(config['reference']),
               "-z",
               "-c",
               "1",
               "-S",
               "2",
               "-E",
               "3",
               "-g",
               "4",
               "-B",
               "{}".format(config['vardict']['num_cores']),
               # "-a", the amplicon flag seems to be creating errors
               # "-F 0", Probably don't need this as duplicates aren't marked and ignoring secondary alignment good
               "-f",
               "{}".format(config['min_alt_af']),
               "-N",
               "{}".format(name),
               "-b",
               "{}".format(input_bam),
               "{}".format(samples[name]['regions'])]

    vardict2vcf = ["{}".format(config['vardict2vcf']['bin']),
                   "-E",
                   "-f",
                   "{}".format(config['min_alt_af']),
                   "-N",
                   "{}".format(name)]

    vcfsort = ["{}".format(config['vcftools_sort']['bin']),
               "-c"]

    command = ("{vardict} | {strandbias} | {vardict2vcf} | "
               "{sort} > {vcf}".format(vardict=" ".join(vardict), strandbias=config['vardict_strandbias']['bin'],
                                       vardict2vcf=" ".join(vardict2vcf), sort=" ".join(vcfsort), vcf=vardict_vcf))

    job.fileStore.logToMaster("VarDict Command: {}\n".format(command))
    run_and_log_command(command, logfile)

    return vardict_vcf


def run_pindel(job, config, name, input_bam):
    """Run Pindel caller for InDel Detection
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str..
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    pindel_config = "{}.pindel_config.txt".format(name)
    output_dir = "{}_pindel".format(name)
    output_vcf = "{}.pindel.vcf".format(name)

    logfile = "{}.pindel.log".format(name)
    vcf_logfile = "{}.pindel2vcf.log".format(name)

    with open(pindel_config, 'w') as bam_config:
        bam_config.write("%s %s %s\n" % (input_bam, config['insert_size'], name))

    command = ("{}".format(config['pindel']['bin']),
               "-f",
               "{}".format(config['reference']),
               "-c",
               "ALL",
               "-w",
               "{}".format(config['pindel']['window']),
               "-E",
               "{}".format(config['pindel']['sensitivity']),
               "-T",
               "{}".format(config['pindel']['num_cores']),
               "-o",
               "{}".format(output_dir),
               "-i",
               "{}".format(pindel_config))

    pindel2vcf_command = ("{}".format(config['pindel2vcf']['bin']),
                          "-r",
                          "{}".format(config['reference']),
                          "-R",
                          "{}".format(config['snpeff']['reference']),
                          "-d",
                          "{}".format(config['snpeff']['reference']),
                          "-he",
                          "0.01",
                          "-G",
                          "-P",
                          "{}".format(output_dir),
                          "-v",
                          "{}".format(output_vcf))

    job.fileStore.logToMaster("Pindel Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("Pindel2vcf Command: {}\n".format(pindel2vcf_command))
    run_and_log_command(" ".join(pindel2vcf_command), vcf_logfile)

    return output_vcf


def vt_normalization(job, config, sample, caller, input_vcf):
    """Decompose and left normalize variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param sample: caller name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.{}.normalized.vcf".format(sample, caller)
    logfile = "{}.{}.vt_normalization.log".format(sample, caller)

    normalization = ["zless",
                     "{}".format(input_vcf),
                     "|",
                     "sed",
                     "'s/ID=AD,Number=./ID=AD,Number=R/'",
                     "|",
                     "{}".format(config['vt']['bin']),
                     "decompose",
                     "-s",
                     "-",
                     "|",
                     "{}".format(config['vt']['bin']),
                     "normalize",
                     "-r",
                     "{}".format(config['reference']),
                     "-",
                     ">",
                     "{}".format(output_vcf)]

    job.fileStore.logToMaster("VT Command: {}\n".format(normalization))
    run_and_log_command(" ".join(normalization), logfile)

    return output_vcf


def merge_variant_calls(job, config, sample, callers, vcf_files):
    """Merge variant calls from multiple variant callers
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param callers: Comma-separated list of VCF callers to tag the ensemble output. Must be in same order as vcf_files.
    :type sample: str.
    :param vcf_files: List of input vcf files for merging.
    :type vcf_files: list.
    :returns:  str -- The output vcf file name.
    """

    merged_vcf = "{}.merged.vcf.gz".format(sample)
    uncompressed_vcf = "{}.merged.vcf".format(sample)
    sorted_vcf = "{}.merged.sorted.vcf".format(sample)

    logfile1 = "{}.merging.log".format(sample)
    logfile2 = "{}.uncompress-merging.log".format(sample)
    logfile3 = "{}.merged_sort.log".format(sample)

    vcf_files_string = " ".join(vcf_files)

    command = ["{}".format(config['ensemble']['bin']),
               "ensemble",
               "-c",
               "{}".format(config['ensemble']['num_cores']),
               "--numpass",
               "1",
               "--names",
               "{}".format(callers),
               "{}".format(merged_vcf),
               "{}".format(config['reference']),
               "{}".format(vcf_files_string)]

    command2 = ["bgzip",
                "-cd",
                "{}".format(merged_vcf),
                ">",
                "{}".format(uncompressed_vcf)]

    command3 = ["{}".format(config['picard']['bin']),
                "SortVcf",
                "SEQUENCE_DICTIONARY={}".format(config['dict']),
                "OUTPUT={}".format(sorted_vcf),
                "INPUT={}".format(uncompressed_vcf)]

    sys.stderr.write("Running commands: \n")
    sys.stderr.write("bcbio-variation-recall Command: {}\n".format(command))
    sys.stderr.write("Uncompression Command: {}\n".format(command2))
    sys.stderr.write("Sort Command: {}\n".format(command3))

    job.fileStore.logToMaster("bcbio-variation-recall Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile1)

    job.fileStore.logToMaster("Uncompression Command: {}\n".format(command2))
    run_and_log_command(" ".join(command2), logfile2)

    job.fileStore.logToMaster("Sort Command: {}\n".format(command3))
    run_and_log_command(" ".join(command3), logfile3)

    # The Index file created by Picard often causes problems with the GATK
    index_file = "{}.idx".format(sorted_vcf)
    os.remove(index_file)

    return sorted_vcf


def annotate_vcf(job, config, name, input_vcf, input_bam):
    """Run GATK's VariantAnnotation on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.annotated.vcf".format(name)
    annotation_logfile = "{}.variantannotation.log".format(name)

    annotation_command = ["{}".format(config['gatk-annotate']['bin']),
                          "-T",
                          "VariantAnnotator",
                          "-R",
                          "{}".format(config['reference']),
                          "-nt",
                          "{}".format(config['gatk-annotate']['num_cores']),
                          "--group",
                          "StandardAnnotation",
                          "--dbsnp",
                          "{}".format(config['dbsnp']),
                          "-I",
                          "{}".format(input_bam),
                          "--variant",
                          "{}".format(input_vcf),
                          "-L",
                          "{}".format(input_vcf),
                          "-o",
                          "{}".format(output_vcf)]

    job.fileStore.logToMaster("GATK VariantAnnotator Command: {}\n".format(annotation_command))
    run_and_log_command(" ".join(annotation_command), annotation_logfile)

    return output_vcf


def filter_variants(job, config, name, input_vcf):
    """Run GATK's VariantFilter on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.filtered.vcf".format(name)
    filter_log = "{}.variantfiltration.log".format(name)

    filter_command = ["{}".format(config['gatk-filter']['bin']),
                      "-T",
                      "VariantFiltration",
                      "-R",
                      "{}".format(config['reference']),
                      "--filterExpression",
                      "'MQ0 > {}'".format(config['mq0_threshold']),
                      "--filterName",
                      "'HighMQ0'",
                      "--filterExpression",
                      "'DP < {}'".format(config['coverage_threshold']),
                      "--filterName",
                      "'LowDepth'",
                      "--filterExpression",
                      "'QUAL < {}'".format(config['var_qual_threshold']),
                      "--filterName",
                      "'LowQual'",
                      "--filterExpression",
                      "'MQ < {}'".format(config['map_qual_threshold']),
                      "--filterName",
                      "'LowMappingQual'",
                      "--variant",
                      "{}".format(input_vcf),
                      "-o",
                      "{}".format(output_vcf)]

    job.fileStore.logToMaster("GATK VariantFiltration Command: {}\n".format(filter_command))
    run_and_log_command(" ".join(filter_command), filter_log)

    return output_vcf


def snpeff(job, config, name, input_vcf):
    """Annotate the specified VCF using snpEff
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.snpEff.{}.vcf".format(name, config['snpeff']['reference'])
    logfile = "{}.snpeff.log".format(name)

    snpeff_command = ["{}".format(config['snpeff']['bin']),
                      "-Xmx{}g".format(config['snpeff']['max_mem']),
                      "-onlyTr {}".format(config['transcripts']),
                      "-v",
                      "{}".format(config['snpeff']['reference']),
                      "{}".format(input_vcf),
                      ">"
                      "{}".format(output_vcf)]

    job.fileStore.logToMaster("snpEff Command: {}\n".format(snpeff_command))
    run_and_log_command(" ".join(snpeff_command), logfile)

    return output_vcf


def vcfanno(job, config, name, samples, input_vcf):
    """Take the specified VCF and use vcfanno to add additional annotations
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.vcfanno.snpEff.{}.vcf".format(name, config['snpeff']['reference'])
    logfile = "{}.vcfanno.log".format(name)

    command = ["{}".format(config['vcfanno']['bin']),
               "-p",
               "{}".format(config['vcfanno']['num_cores']),
               "--lua",
               "{}".format(config['vcfanno']['lua']),
               "{}".format(samples[name]['vcfanno_config']),
               "{}".format(input_vcf),
               ">",
               "{}".format(output_vcf)]

    job.fileStore.logToMaster("VCFAnno Command: {}\n".format(command))
    run_and_log_command(" ".join(command), logfile)

    return output_vcf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-a', '--address', help="IP Address for Cassandra connection", default='127.0.0.1')
    parser.add_argument('-u', '--username', help='Cassandra username for login', default=None)
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    cwd = os.getcwd()
    sys.stdout.write("Setting up analysis directory\n")

    # if not os.path.exists("Logs"):
    #     os.makedirs("Logs")
    # if not os.path.exists("FinalVCFs"):
    #     os.makedirs("FinalVCFs")
    # if not os.path.exists("FinalBAMs"):
    #     os.makedirs("FinalBAMs")
    # if not os.path.exists("Intermediates"):
    #     os.makedirs("Intermediates")
    # if not os.path.exists("Coverage"):
    #     os.makedirs("Coverage")
    # if not os.path.exists("Reports"):
    #     os.makedirs("Reports")

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    if args.username:
        password = getpass.getpass()
        auth_provider = PlainTextAuthProvider(username=args.username, password=password)
    else:
        auth_provider = None

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(spawn_batch_jobs, cores=1)

    fastqc_job = Job.wrapJobFn(run_fastqc, config, samples)

    # Per sample jobs
    for sample in samples:
        # Alignment and Refinement Stages
        align_job = Job.wrapJobFn(run_bwa_mem, config, sample, samples,
                                  cores=int(config['bwa']['num_cores']),
                                  memory="{}G".format(config['bwa']['max_mem']))

        add_job = Job.wrapJobFn(add_or_replace_readgroups, config, sample, align_job.rv(),
                                cores=1,
                                memory="{}G".format(config['picard-add']['max_mem']))

        creator_job = Job.wrapJobFn(realign_target_creator, config, sample, add_job.rv(),
                                    cores=int(config['gatk-realign']['num_cores']),
                                    memory="{}G".format(config['gatk-realign']['max_mem']))

        realign_job = Job.wrapJobFn(realign_indels, config, sample, add_job.rv(), creator_job.rv(),
                                    cores=1,
                                    memory="{}G".format(config['gatk-realign']['max_mem']))

        recal_job = Job.wrapJobFn(recalibrator, config, sample, realign_job.rv(),
                                  cores=int(config['gatk-recal']['num_cores']),
                                  memory="{}G".format(config['gatk-recal']['max_mem']))

        coverage_job = Job.wrapJobFn(sambamba_region_coverage, config, sample, samples,
                                     "{}.recalibrated.sorted.bam".format(sample),
                                     cores=int(config['gatk']['num_cores']),
                                     memory="{}G".format(config['gatk']['max_mem']))

        # Variant Calling
        spawn_variant_job = Job.wrapJobFn(spawn_batch_jobs)

        freebayes_job = Job.wrapJobFn(freebayes_single, config, sample,
                                      "{}.recalibrated.sorted.bam".format(sample),
                                      cores=1,
                                      memory="{}G".format(config['freebayes']['max_mem']))

        mutect_job = Job.wrapJobFn(mutect_single, config, sample, samples,
                                   "{}.recalibrated.sorted.bam".format(sample),
                                   cores=1,
                                   memory="{}G".format(config['mutect']['max_mem']))

        vardict_job = Job.wrapJobFn(vardict_single, config, sample, samples,
                                    "{}.recalibrated.sorted.bam".format(sample),
                                    cores=int(config['vardict']['num_cores']),
                                    memory="{}G".format(config['vardict']['max_mem']))

        scalpel_job = Job.wrapJobFn(scalpel_single, config, sample, samples,
                                    "{}.recalibrated.sorted.bam".format(sample),
                                    cores=int(config['scalpel']['num_cores']),
                                    memory="{}G".format(config['scalpel']['max_mem']))

        platypus_job = Job.wrapJobFn(platypus_single, config, sample, samples,
                                     "{}.recalibrated.sorted.bam".format(sample),
                                     cores=int(config['platypus']['num_cores']),
                                     memory="{}G".format(config['platypus']['max_mem']))

        pindel_job = Job.wrapJobFn(run_pindel, config, sample,
                                   "{}.recalibrated.sorted.bam".format(sample),
                                   cores=int(config['pindel']['num_cores']),
                                   memory="{}G".format(config['pindel']['max_mem']))

        # Need to filter for on target only results somewhere as well
        spawn_normalization_job = Job.wrapJobFn(spawn_batch_jobs)

        normalization_job1 = Job.wrapJobFn(vt_normalization, config, sample, "freebayes",
                                           "{}.freebayes.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job2 = Job.wrapJobFn(vt_normalization, config, sample, "mutect",
                                           "{}.mutect.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job3 = Job.wrapJobFn(vt_normalization, config, sample, "vardict",
                                           "{}.vardict.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job4 = Job.wrapJobFn(vt_normalization, config, sample, "scalpel",
                                           "{}.scalpel.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job5 = Job.wrapJobFn(vt_normalization, config, sample, "platypus",
                                           "{}.platypus.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job6 = Job.wrapJobFn(vt_normalization, config, sample, "pindel",
                                           "{}.pindel.vcf".format(sample),
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        callers = "freebayes,mutect,vardict,scalpel,platypus,pindel"

        merge_job = Job.wrapJobFn(merge_variant_calls, config, sample, callers, (normalization_job1.rv(),
                                                                                           normalization_job2.rv(),
                                                                                           normalization_job3.rv(),
                                                                                           normalization_job4.rv(),
                                                                                           normalization_job5.rv(),
                                                                                           normalization_job6.rv()))

        gatk_annotate_job = Job.wrapJobFn(annotate_vcf, config, sample, merge_job.rv(),
                                          "{}.recalibrated.sorted.bam".format(sample),
                                          cores=int(config['gatk-annotate']['num_cores']),
                                          memory="{}G".format(config['gatk-annotate']['max_mem']))

        gatk_filter_job = Job.wrapJobFn(filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1,
                                        memory="{}G".format(config['gatk-filter']['max_mem']))

        snpeff_job = Job.wrapJobFn(snpeff, config, sample, "{}.filtered.vcf".format(sample),
                                   cores=int(config['snpeff']['num_cores']),
                                   memory="{}G".format(config['snpeff']['max_mem']))

        vcfanno_job = Job.wrapJobFn(vcfanno, config, sample, samples,
                                    "{}.snpEff.{}.vcf".format(sample, config['snpeff']['reference']),
                                    cores=int(config['vcfanno']['num_cores']),
                                    memory="{}G".format(config['vcfanno']['max_mem']))

        # Create workflow from created jobs
        root_job.addChild(align_job)
        align_job.addChild(add_job)
        add_job.addChild(creator_job)
        creator_job.addChild(realign_job)
        realign_job.addChild(recal_job)

        recal_job.addChild(spawn_variant_job)

        spawn_variant_job.addChild(coverage_job)
        spawn_variant_job.addChild(freebayes_job)
        spawn_variant_job.addChild(mutect_job)
        spawn_variant_job.addChild(vardict_job)
        spawn_variant_job.addChild(scalpel_job)
        spawn_variant_job.addChild(platypus_job)
        spawn_variant_job.addChild(pindel_job)

        spawn_variant_job.addFollowOn(spawn_normalization_job)

        spawn_normalization_job.addChild(normalization_job1)
        spawn_normalization_job.addChild(normalization_job2)
        spawn_normalization_job.addChild(normalization_job3)
        spawn_normalization_job.addChild(normalization_job4)
        spawn_normalization_job.addChild(normalization_job5)
        spawn_normalization_job.addChild(normalization_job6)

        spawn_normalization_job.addFollowOn(merge_job)

        merge_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(snpeff_job)
        snpeff_job.addChild(vcfanno_job)

    root_job.addFollowOn(fastqc_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
