% function cvsets_nestfamily()
% function to set up cross-validation to respect family membership

clear all; close all

% data dir
%data_dir = 'input_data/';
%data_dir = 'input_data_noBlueBrains/';
data_dir = 'input_data_noBlueBrains_parcellationfree/'; %slightly fewer subs to ensure all had 1 border and 1 ectopic

% variables
cv_num = 10; % number of crossvalidation subsets

% load subject data
load([data_dir 'SubjectInfo_forCaterina.mat']);
subject_info_labels{6} = 'sex'; %quick fix for annoying variable name

% put into a datatable
subInfoTable = cell2table(subject_info,'VariableNames',subject_info_labels);

% figure out "families" - i.e., groups that share EITHER a mom or a dad ID
famnum = 1;
familyID = nan*ones(size(subInfoTable,1),1);
for s = 1:size(subInfoTable,1)
    thisMom = subInfoTable.momID(s);
    thisDad = subInfoTable.dadID(s);

    if isnan(familyID(s)) % give them a new family id, along with everyone related to them
        familyID(subInfoTable.momID == thisMom) = famnum;
        familyID(subInfoTable.dadID == thisDad) = famnum;
        famnum = famnum+1;
    else
        % double check that everyone related to them has the same family ID
        anymissing1 = sum(familyID(subInfoTable.momID == thisMom) ~= familyID(s));
        anymissing2 = sum(familyID(subInfoTable.dadID == thisDad) ~= familyID(s));
        if (anymissing1 + anymissing2) > 0
            disp('WARNING: missing family members');
        end
    end
end

% CREATE CV SUBSETS 

% first, find all families
unique_fams = unique(familyID);
num_unique_fams = length(unique_fams);
% WAS going to go in order of size to limit remainders
% BUT realized that would make initial cv sets have all the biggest
% families - not ideal.
% for f = 1:num_unique_fams 
%     unique_fam_size(f) = sum(strcmp(unique_fams{f},subInfoTable.familyID));
% end
% [sorted_fam_size sorted_fam_idx] = sort(unique_fam_size,'descend');
% unique_fams_sorted = unique_fams(sorted_fam_idx);

% permute unique_fams order once so not in numerical order
rng(987); %rng(3); %make this predictable for now
unique_fams_sorted = unique_fams(randperm(length(unique_fams)));

% now go through each cv and add family groups to it until its full
fami = 1; % initial index
num_sub_percv = floor(size(subInfoTable,1)/cv_num); % left overs will be added at end
for i = 1:cv_num
    
    % initialize the array
    cv_set{i} = [];
    cv_setIDs{i} = [];
    cv_not_full = 1;
    
    while cv_not_full
        
        % loop through families and select subs
        theseSubs = find(unique_fams_sorted(fami)==familyID);

        % updated cv length
        cv_length_new = length(cv_set{i}) + length(theseSubs);
        
        % update cv set with this family if it can fit
        if cv_length_new <= num_sub_percv
            cv_set{i}(end+1:end+length(theseSubs)) = subInfoTable.subjectID(theseSubs);
            cv_setIDs{i}(end+1:end+length(theseSubs)) = theseSubs;
            fami = fami+1; %update family number
        else
            cv_not_full = 0; %break out of loop
            cv_set_size(i) = length(cv_set{i});
            
        end
    end
end

% check for leftovers and apportion as possible
% note, right now only missing one is the nan group of 10 which is hard to
% place at the end, otherwise groups reasonably close in size. Leave as is?
% allow slight differences in size, lumping final remainder group into
% smallest cv?
while num_unique_fams >= fami
    disp('MISSING SOME SUBS IN final sets - placing in smallest cv group');
    [val, smallest_cv] = min(cv_set_size);
    
    % update subs
    theseSubs = find(unique_fams_sorted(fami)==familyID);
    cv_set{smallest_cv}(end+1:end+length(theseSubs)) = subInfoTable.subjectID(theseSubs);
    cv_setIDs{smallest_cv}(end+1:end+length(theseSubs)) = theseSubs;
    
    %update other values
    cv_set_size(smallest_cv) = length(cv_set{smallest_cv});
    
    fami = fami+1;
end

% quick check on how sets are balanced in terms of age and sex
for i = 1:cv_num
    cv_avg_age(i,1) = mean(subInfoTable.age(cv_setIDs{i}));
    cv_avg_age(i,2) = std(subInfoTable.age(cv_setIDs{i}))/sqrt(length(cv_setIDs{i}));
    cv_avg_sex(i,1) = mean(subInfoTable.sex(cv_setIDs{i}));
    cv_avg_sex(i,2) = std(subInfoTable.sex(cv_setIDs{i}))/sqrt(length(cv_setIDs{i}));
end

% write it out in a different format that will be easier for SVR edits
% (i.e., to remove beh subs missing data)
subs_byCVnum = nan*ones(size(subInfoTable,1),1);
for i = 1:cv_num
    subs_byCVnum(cv_setIDs{i}) = i;
end

save([data_dir 'HCP_CVsets_nestedfamily.mat'],'subs_byCVnum','cv_setIDs');