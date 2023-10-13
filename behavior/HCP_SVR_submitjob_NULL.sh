#!/bin/bash
#SBATCH --account=p31149 ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal ### PARTITION (buyin, short, normal, w10001, etc)
#SBATCH --array=0-10 ##57 ## number of jobs to run "in parallel"
#SBATCH --nodes=1 ## Never need to change this
#SBATCH --ntasks-per-node=1 ## Never need to change this
#SBATCH --time=47:59:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=40G ## how much RAM you need per computer (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name="HCP_SVR_job_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=SVR_YeoLocT1000_jobNULL.%A_%a.out ## use the jobid (A) and the specific job index (a) to name your log file

# read in each row of the input_arguments text file into an array called input_args
IFS=$'\n' read -d '' -r -a input_arguments < SVR_nodes_inp_YeoBeh_location_TEMP.txt

# Use the SLURM_ARRAY_TASK_ID variable to select the correct index from input_arguments and then split the string by whatever delimiter you chose (in our case each input argument is split by a space)
IFS=' ' read -r -a input_args <<< "${input_arguments[$SLURM_ARRAY_TASK_ID]}"

# Pass the input arguments associated with this SLURM_ARRAY_TASK_ID to your function.
## job commands; job_array is the MATLAB .m file, specified without the .m extension
module load matlab/r2020b

# matlab -singleCompThread -batch "input_arg1='YeoBeh';input_arg2='network';input_arg3='border';input_arg4='2';HCP_SVR_analysis_nodes"
matlab -singleCompThread -batch "input_arg1='${input_args[0]}';input_arg2='${input_args[1]}';input_arg3='${input_args[2]}';input_arg4='${input_args[3]}';HCP_SVR_analysis_nodes_null"
#matlab -singleCompThread -batch "HCP_SVR_analysis_nodes_null(${input_args[0]},${input_args[1]},${input_args[2]},${input_args[3]})"
