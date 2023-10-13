function [behScores, behLabels] = load_behavioral_data(input_dir,beh_type)

switch beh_type
    case 'YeoBeh'
        load([input_dir 'YeoBeh.mat']); % loads Yeo behavioral variables
        behScores = YeoBehScores;
        %behLabels = YeoBehLabels;
        behLabels = YeoBehLabelsLong;
    case 'SeitzmanBeh'
        load([input_dir 'SeitzmanBeh.mat']);
        behScores = SeitzmanBehScores;
        behLabels = SeitzmanBehLabels;
        %behLabelsLong = SeitzmanBehLabelsLong;
    otherwise
        error('Do not recognize behavioral label type');
end

end