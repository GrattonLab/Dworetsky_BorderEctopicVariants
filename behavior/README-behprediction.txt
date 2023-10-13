Instructions for running

1. Run SVR prediction on actual data

Runs on cluster using:
> sbatch HCP_SVR_submitjob.sh

This calls on
> HCP_SVR_analysis_nodes
with inputs specified in:
> SVR_nodes_inp_YeoBeh_{network/location}_{border/ectopic}.txt (set this correctly in the script before batch command call!)

Which is one per line of:
YeoBeh (other option was SeitzmanBeh which gave similar results but many fewer variables, so harder to tell relationships)
network vs. location in turn (prediction based on network affiliation or location of variants)
border vs. ectopic in turn (prediction based on border vs. ectopic variants)
Each run separately in parallel on each behavioral variable (one row on txt file)

Requires the following other functions:
load_subject_data.m : load subject demographic info
load_fMRI_data.m : loads fMRI data for each type of analysis and reformats as needed
load_behavioral_data.m : loads subject behavioral data of appropriate type, reformats as needed
tenFOLD_svm_Dworetsky.m : wrapper script for SVM from A. Nielsen, setting up cross-validation and saving out relevatn variables
svm_scripts_v2021_Dworetsky.m: script that actually does the SVM fit and test on a particular split of the data

Note some paths set at the top of the matlab script which will need to be changed if run elsewhere
You should be able to just edit "top_dir", but make sure output_dirs are already created with correct names

Some notes:
- if 'beh_regress' is set to true (as in paper), regression of the following nuisance behavioral variables is performed: age, sex, FD, DVARS, BMI (done within each CV and applied to left out data)
- a pre-specified set of CV groups was used to ensure that familial structure was respected in CV (created through cvsets_nestfamily_updated.m function)
- features representing correlations are r-to-z transformed
- location based version of prediction takes a fair amount of memory (and more time) to run

To run without cluster call for testing:
% set workspace variables and call
>> input_arg1 = 'YeoBeh';
>> input_arg2 = 'network';
>> input_arg3 = 'border';
>> input_arg4 = '1';
>> HCP_SVR_analysis_nodes

Notes on running through sbatch for these jobs on our machine:
- network versions are fast and don't require much memory (<20G, short queue, done in a minute)
- location versions are slower and require more memory  (<40G, short queue, done in 10 minutes or so)
- nulls take a long time (XXX 1000x each of the versions above)

-----------------------------------------------------------
2. Run SVR prediction on null data
> HCP_SVR_submitjob_NULL.sh

same as above, but calls on:
> HCP_SVR_analysis_nodes_null.m
which specifies nperms (500)
and permutes behavioral variable order before each permutation

takes a long time to run

------------------------------------------------------------
3. Plot final data, output analyses:
> HCP_SVR_nodes_finalplotting.m 

Both #1 and #2 need to be run for this to work
(although can proceed without #2 and just omit plotting of permutation lines)

depends on:
save_fig : Grattonlab function to save out plots easily; dependent on operating system, will likely need to be altered for specific users