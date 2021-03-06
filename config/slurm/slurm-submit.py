#!/usr/bin/env python
import argparse
import subprocess
import sys

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "positional", action="append",
    nargs="?", metavar="POS",
    help="additional arguments not in slurm parser group to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "-A", "--account", help="charge job to specified account")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task")
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-e", "--error", help="file for batch script's standard error")
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks", help="number of tasks to run")
slurm_parser.add_argument(
    "-N", "--nodes", help="number of nodes on which to run (N = min[-max])")
slurm_parser.add_argument(
    "-o", "--output", help="file for batch script's standard output")
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "--gres", help="generic resources")
slurm_parser.add_argument(
    "--gres-flags", help="flags related to GRES")
slurm_parser.add_argument(
    "--wrap", help="wrap command string in a sh script and submit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")
slurm_parser.add_argument(
    "--mem", help="minimum amount of real memory")

args = parser.parse_args()

if args.help:
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

extras = ""
if args.positional:
    for m in args.positional:
        if m is not None:
            extras = extras + " " + m

arg_dict = dict(args.__dict__)

# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    arg_dict["time"] = "{}:00:00".format(resources.get("time_hr", arg_dict["time"]))
    arg_dict["mem"] = '{}g'.format(resources.get("mem_gb", arg_dict["mem"]))

# Threads
if "threads" in job_properties:
    arg_dict["ntasks"] = job_properties["threads"]

# Set default partition
if arg_dict["partition"] is None:
    arg_dict["partition"] = "norm"

# Add quick if walltime is less than 4hrs
if int(arg_dict["time"].split(':')[0]) <= 4:
    arg_dict["partition"] = "quick" + "," + arg_dict["partition"]

opt_keys = ["array", "account", "begin", "cpus-per-task", "mem",
            "partition", "depedency", "workdir", "error", "job-name", "mail-type",
            "mail-user", "ntasks", "nodes", "output", "quiet", "time",
            "wrap", "constraint", "gres", "gres-flags"]
opts = ""
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opts += " --{} \"{}\" ".format(k, v)

if arg_dict["wrap"] is not None:
    cmd = "sbatch {opts}".format(opts=opts)
else:
    cmd = "sbatch {opts} {extras}".format(opts=opts, extras=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
print(res)
