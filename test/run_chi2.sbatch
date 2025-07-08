#!/bin/bash
#SBATCH --job-name=Chi2
#SBATCH --output=results/out_%j.log
#SBATCH --error=results/err_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

source /store/sw/anaconda3/etc/profile.d/conda.sh
conda activate ds4hep

echo "---------- Starting Hypothesis tests ----------"
python ./run_chi2.py -s ${V_SIN14} -m ${V_DM41} -n ${V_NSIG} -o ${V_FNAME}
echo "-----------------------------------------------"

