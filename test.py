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
        create_job = Job.wrapJobFn(create_file, fname)

        root_job.addChild(create_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
