#!/bin/bash

#SBATCH --job-name=ArcSimPerm_JB22_Bcells
#SBATCH --output=slurm/ArcSimPerm_JB22_Bcells/%j.out
#SBATCH --error=slurm/ArcSimPerm_JB22_Bcells/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=72:00:00

# Nextflow
pixi run nextflow run nextflow/pipelines/ArcSimPerm_JB22_Bcells.nf \
    -c nextflow/configs/ArcSimPerm_JB22_Bcells.config \
    -work-dir /mnt/ceph/groups_hdd/SCCGroup/computational_systems_biology_group/dkatzlberger/ArcSimulationPermutation/work/ArcSimPerm_JB22_Bcells \
    -resume