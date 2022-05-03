#!/usr/bin/env python3
#makePGLSbatch.py [params] > foo.sb

import argparse


parser = argparse.ArgumentParser(description="Generate SBATCH file for PGLS", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-p", required=False, help="String of partitions to use", default="pfen_bigmem,pfen1,pfen2,pfen3,pool1")
parser.add_argument("-j", help="job name", default="pgls")
parser.add_argument("-l", help="log file directory", default="/home/dschaffe/enhancer-clustering/scripts/log/")
parser.add_argument("-r", default=4000, help="Memory to request in MB", type=int)
parser.add_argument("-t", default="0-8", help="time limit per job (Slurm format)")
parser.add_argument("-n", help="number of array jobs", type=int, required=True)
parser.add_argument("-z", help="Slurm-format list of array jobs (default 1-n)")
parser.add_argument("-s", help="number of permulations (0=one iteration w/o permulation)", default=0, type=int)
parser.add_argument("-d", help="Base data directory", default="/home/dschaffe/enhancer-clustering/data/")
parser.add_argument("-y", help="Tree file (relative to data)", default="200mammals_x_lessgc40_75consensus.tree")
parser.add_argument("-m", help="Matrix of activity predictions (relative to data)", required=True)
parser.add_argument("-c", help="List of species (column names) for matrix", default="BoreoeutheriaSpeciesNamesNew.txt")
parser.add_argument("-a", help="Zoonomia-format file of annotations (relative to data)", required=True)
parser.add_argument("-b", help="Column name of phenotype", required=True)
parser.add_argument("-o", help="Output file template (relative to data)", required=True)
parser.add_argument("-q", help="Method to use (phyloglm, phylolm, etc.)", required=True)
parser.add_argument("-e", help="Path to enhancer_[method].r", required=False, default="/home/dschaffe/enhancer-clustering/scripts")


args = parser.parse_args()

print("#!/bin/bash")
print("#SBATCH --partition=" + args.p)
print("#SBATCH --job-name=" + args.j)
print("#SBATCH --cpus-per-task=1") 
print("#SBATCH --error=" + args.l + "/" + args.j + "_%A_%a.err.txt")
print("#SBATCH --output=" + args.l + "/" + args.j + "_%A_%a.out.txt")
print("#SBATCH --mem=" + str(args.r) + "M") 
print("#SBATCH --time=" + args.t)
arr = args.z if args.z else "1-" + str(args.n)
print("#SBATCH --array=" + arr)
print()
print("seed=`od --read-bytes=4 --address-radix=n --format=u4 /dev/random | awk '$1>=2^31{$1-=2^32}1'`")
print("data=\"" + args.d + "\"")
print("tree=\"" + args.d + "/" + args.y + "\"")
print("preds=\"" + args.d + "/" + args.m + "\"")
print("species=\"" + args.d + "/" + args.c + "\"")
print("traits=\"" + args.d + "/" + args.a + "\"")
print("out=\"" + args.d + "/" + args.o + "\"")

print("Rscript " + args.e + "/enhancer_" + args.q + ".r $tree $preds $species $traits $out ${SLURM_ARRAY_TASK_ID}", int(args.n), int(args.s), \
        "$seed", args.b, args.e)
