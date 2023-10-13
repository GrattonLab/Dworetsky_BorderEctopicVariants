function subInfoTable = load_subject_data(input_dir)

load([input_dir '/SubjectInfo_forCaterina.mat']);
subject_info_labels{6} = 'sex'; %quick fix for annoying variable name
subInfoTable = cell2table(subject_info,'VariableNames',subject_info_labels);

end