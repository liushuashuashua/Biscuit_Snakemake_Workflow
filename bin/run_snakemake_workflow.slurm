#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J SNAKEMASTER
#SBATCH -o logs/workflows/workflow_output-%j.log
#SBATCH --ntasks 1
#SBATCH --time 240:00:00
#SBATCH --mem=8G

mkdir -p logs/workflows

cd ${SLURM_SUBMIT_DIR}

# Add snakemake to PATH here
if [[ `which snakemake 2>&1 /dev/null` ]]; then
    snakemake_module="bbc2/snakemake/snakemake-7.25.0"

    module load $snakemake_module
fi

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --configfile config/config.yaml --dry-run         > logs/workflows/workflow_${TIME}.txt
snakemake --configfile config/config.yaml --dag | dot -Tpng > logs/workflows/workflow_${TIME}.png

# Default to using conda, if using environment modules, then replace --use-conda with --use-envmodules
# Note, this requires downloading mamba (conda install -n base -c conda-forge mamba)
snakemake \
    --printshellcmds \
    --latency-wait 20 \
    --keep-going \
    --use-conda \
    --jobs 20 \
    --configfile config/config.yaml \
    --cluster "mkdir -p logs/{rule}; sbatch \
        --export=ALL \
        --ntasks {threads} \
        --mem={resources.mem_gb}G \
        -t {resources.time} \
        -o logs/{rule}/{rule}-{wildcards}-%j.log"
