This repository contains MATLAB scripts used to conduct border/ectopic variant analyses from Dworetsky et al., 2024, *Two common and distinct forms of variation in human functional brain networks* (Nature Neuroscience). These scripts allow users to classify network variants as either border shifts or ectopic intrusions, and compare their spatial/network distributions, prediction of behavioral phenotypes, activation during tasks, subgroup properties, and genetic similarity.

For original variant creation steps, see scripts in GitHub repo from Seitzman et al., 2019, *Trait-like variants in human functional brain networks* (PNAS): https://github.com/GrattonLab/SeitzmanGratton-2019-PNAS to compute individual-to-group spatial correlation, binarize pre-variant maps, and extract variants' correlation coefficients to network templates.

1. **variant_classification**: Classify variants as border or ectopic using classify_variants_network_dependent.m (parcellation-dependent, distance-to-same-network classification) or classify_variants_parcellation_free.m (parcellation-independent, distance-to-peak-group-similarity classification).

The scripts in remaining folders can be run in any order:
- **spatial_distribution**: perform and visualize results of permutation analysis examining border/ectopic variant differences in spatial distribution across the cortex, plus cluster-correction of permutation analyses.
- **network_distribution**: perform and plot results of permutation analysis examining border/ectopic variant differences in network assignments.
- **twins**: compute and plot Dice coefficient of subject pairs of interest across groups (e.g., MZ twins, DZ twins, siblings, unrelated individuals); permutation analysis of Falconer's formula in MZ and DZ groups.
- **task responses**: organize task contrast map data; compare border and ectopic task activations across contrasts in relation to canonical assigned network and all other networks; plot normalized value of variants' activations shift toward canonical network (for networks of interest).
- **subgroups**: identify/plot properties of border/ectopic subgroups across individuals. Intended to be run after templateMatchingVariants.m from https://github.com/GrattonLab/SeitzmanGratton-2019-PNAS (separately for border and ectopic variants). Also requires Infomap (https://www.mapequation.org) - Run_Infomap_2015 has been previously shared at https://github.com/MidnightScanClub/MSCcodebase/tree/master/Utilities/Infomap_wrapper.
- **behavior**: use variants' network or location info to predict behavioral phenotypes (see README-behprediction.txt for more info).

**Resources**:

Each of these scripts requires supporting scripts for reading and writing CIFTI files (see more at https://github.com/fieldtrip/fieldtrip/). We have included modified versions of these scripts in the *resources* folder. This folder also includes a scripts for computing Dice coefficients and a script for plotting jittered distributions.
