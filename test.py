#!/usr/bin/env python

# Standard packages
import os
import sys
import argparse
import configuration
import subprocess as sub

from toil.job import Job
from toil.job import JobException


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")


def create_file(job, fpath):
    job.fileStore.logToMaster("Creating file: {}\n".format(fpath))
    sub.call(['touch', fpath])
    # job.fileStore.writeGlobalFile(fpath)


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
               "24",
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
                "24",
                "-O",
                "bam",
                "-o",
                "{}".format(output_bam),
                "-T",
                "{}".format(temp),
                "-"]

    command = "{} | {} | {}".format(" ".join(bwa_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("BWA Command: {}\n".format(command))

    with open(logfile, "wb") as err:
        sys.stdout.write("Executing {} and writing to logfile {}\n".format(command, logfile))
        err.write("Command: {}\n".format(command))
        # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        # output = p.communicate()
        # code = p.returncode
        # if code:
        #     raise RuntimeError("An error occurred when executing the commandline: {}. "
        #                        "Please check the logfile {} for details\n".format(command, logfile))

    return output_bam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    cwd = os.getcwd()

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    root_job = Job.wrapJobFn(spawn_batch_jobs, cores=1)

    for sample in samples:
        fname = "{}/{}.txt".format(cwd, sample)
        align_job = Job.wrapJobFn(run_bwa_mem, config, sample, samples, cores=24,
                                  memory="112G".format(config['bwa']['max_mem']))

        root_job.addChild(align_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
