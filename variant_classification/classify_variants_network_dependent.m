function classify_variants_network_dependent(variantsInfoFile, consensusNetsFile, vertexNeighborsFile, distances, outdir)

%
% This function classifies variants as border or ectopic based on their
% minimum (edge-to-edge) distance between variant vertices and same-network
% vertices in the group consensus network map. The codebase for producing
% necessary inputs (variants labeled by ID and by network) can be found at:
% https://github.com/GrattonLab/SeitzmanGratton-2019-PNAS
% This script requires supporting scripts for reading and writing CIFTI
% files, where the modified versions contained in this script can be found
% at the above link.
%
% INPUTS:
% variantsInfoFile: a text file with 3 columns - (1) subject ID; (2) path to 
%   file containing that subject's variants, labeled by a unique ID starting 
%   at 1;(3) path to file containing that subject's variants labeled with
%   numbers corresponding to assigned networks
% consensusNetsFile: path to a CIFTI with map of group consensus networks
% vertexNeighborsFile: path to .mat file containing each vertex's neighboring
%   vertices
% distances: a node-to-node distance matrix across all vertices
% outdir: output directory for border/ectopic classified variants
%
% OUTPUTS:
% a CIFTI file for each subject, where border variant locations are 
%   labeled with a value of 1 and ectopic variant locations are labeled 
%   with a value of 2
%


%%% change if desired %%%
ectopicDistance = 3.5;  %in millimeters

% load, read, and set up variables
[subjects, variantsIDs, variantsAssigned] = textread(variantsInfoFile, '%s%s%s');
consen = ft_read_cifti_mod(consensusNetsFile); consen = consen.data;
load(vertexNeighborsFile);

% loop through subjects
for sub = 1:length(subjects)
    
    subID = subjects{sub};
    disp(['Classifying variants for subject ' subID ' (' num2str(sub) '/' num2str(length(subjects)) ')'])
    
    % load variants labeled by ID and by network assignment, set up border/ectopic data file
    varsUnique = ft_read_cifti_mod(variantsIDs{sub});
    varsUnique = varsUnique.data;
    varsAssign = ft_read_cifti_mod(variantsAssigned{sub});
    b1e2 = varsAssign; b1e2.data = zeros(size(varsAssign.data,1), 1);
    varsAssign = varsAssign.data;

    subVarIDs = unique(varsUnique); subVarIDs(subVarIDs==0)=[];
    
    for var = 1:length(subVarIDs)
        varInds = find(varsUnique == subVarIDs(var));
        varNet = varsAssign(varInds(1));
        
        % check whether each variant vertex shares neighbors with same
        % network in group map
        consInds = find(consen == varNet);
        neighboringInds = 0;
        for ind = 1:length(varInds)
            nbs = neighbors(varInds(ind),2:end); nbs(isnan(nbs)) = [];
            for neigh = 1:length(nbs)
                if ismember(nbs(neigh),consInds)
                    neighboringInds = neighboringInds + 1;
                end
            end    
        end
        
        % if any overlap, variant is a border variant
        if neighboringInds > 0
            b1e2.data(varInds) = 1;
        
        % if no overlap, check distance from nearest same-network vertices
        % in the group map; if < [ectopicDistance] mm, variant is a border
        % variant; otherwise, variant is ectopic
        else
            dists = [];
            for i = 1:length(varInds)
                dists = [dists; distances(varInds(i), consInds(:))];
            end
            mindist = min(min(dists));
            if mindist <= ectopicDistance
                b1e2.data(varInds) = 1;
            else
                b1e2.data(varInds) = 2;
            end
        end
        
        clear varInds varNet consInds
    end
    
    ft_write_cifti_mod([outdir '/' subID '_border1ectopic2.dtseries.nii'], b1e2)
    
    clear varsUnique varsAssign numVars b1e2
end

end
