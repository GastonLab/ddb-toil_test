#!/usr/bin/env python

# Standard packages
import os
import sys
import argparse
import configuration
import subprocess

from toil.job import Job
from toil.job import JobException


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")


def create_file(job, fpath):
    job.fileStore.logToMaster("Creating file: {}\n".format(fpath))
    subprocess.call(['touch', fpath])
    # job.fileStore.writeGlobalFile(fpath)


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
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            raise RuntimeError("An error occurred when executing the commandline: {}. "
                               "Please check the logfile {} for details\n".format(command, logfile))


def run_bwa_mem(job, workdir, config, sample, samples):
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

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(sample))

    output_bam = "{}.bwa.sorted.bam".format(sample)
    temp = "{}.bwa.sort.temp".format(sample)
    logfile = "{}.bwa-align.log".format(sample)

    cmd = ["{}".format(config['bwa']['bin']),
           "mem",
           "-t",
           "24",
           "-M",
           "-v",
           "2",
           "{}".format(config['reference']),
           "{}".format(samples[sample]['fastq1']),
           "{}".format(samples[sample]['fastq2']),
           ">",
           "{}.aligned.sam".format(sample)]

    command = " ".join(cmd)
    job.fileStore.logToMaster("BWA Command: {}\n".format(command))
    subprocess.check_call(command, shell=True)
    return job.fileStore.writeGlobalFile(os.path.join(workdir, "{}.aligned.sam".format(sample)))


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
        align_job = Job.wrapJobFn(run_bwa_mem, cwd, config, sample, samples,
                                  cores=24,
                                  memory="112G".format(config['bwa']['max_mem']))

        root_job.addChild(align_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
