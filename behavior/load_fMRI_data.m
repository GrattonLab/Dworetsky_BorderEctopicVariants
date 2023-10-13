function thisVariantInput = load_fMRI_data(input_dir, variant_form, variant_type)

switch variant_form
    case 'network'
        load([input_dir 'networkInfo.mat']); % loads info about network corrs 
        switch variant_type
            case 'border'
               borderNetworkCorrs = borderNetworkCorrs'; % go to 11x814
               thisVariantInput = reshape(borderNetworkCorrs,[size(borderNetworkCorrs,1) 1 size(borderNetworkCorrs,2)]); %11 x 1 x 814
            case 'ectopic'
               ectopicNetworkCorrs = ectopicNetworkCorrs'; % go to 11x814
               thisVariantInput = reshape(ectopicNetworkCorrs,[size(ectopicNetworkCorrs,1) 1 size(ectopicNetworkCorrs,2)]); %11 x 1 x 814
            case 'all'
                borderNetworkCorrs = borderNetworkCorrs'; % go to 11x814
                thisVariantInput_border = reshape(borderNetworkCorrs,[size(borderNetworkCorrs,1) 1 size(borderNetworkCorrs,2)]); %11 x 1 x 814
                ectopicNetworkCorrs = ectopicNetworkCorrs'; % go to 11x814
                thisVariantInput_ectopic = reshape(ectopicNetworkCorrs,[size(ectopicNetworkCorrs,1) 1 size(ectopicNetworkCorrs,2)]); %11 x 1 x 814
                thisVariantInput = [thisVariantInput_border; thisVariantInput_ectopic];
        end
    case 'location'
        load([input_dir 'locationInfo.mat']); %loads info about variant locations
        switch variant_type
            case 'border'
                thisVariantInput = reshape(allbordervars,[size(allbordervars,1) 1 size(allbordervars,2)]); %starts 59412 x 814, change to 59412 x 1 x 814
            case 'ectopic'
                thisVariantInput = reshape(allectopicvars,[size(allectopicvars,1) 1 size(allectopicvars,2)]); %starts 59412 x 814, change to 59412 x 1 x 814
            case 'all'
                thisVariantInput_border = reshape(allbordervars,[size(allbordervars,1) 1 size(allbordervars,2)]); %starts 59412 x 814, change to 59412 x 1 x 814
                thisVariantInput_ectopic = reshape(allectopicvars,[size(allectopicvars,1) 1 size(allectopicvars,2)]); %starts 59412 x 814, change to 59412 x 1 x 814
                %thisVariantInput = thisVariantInput_border + thisVariantInput_ectopic; % sum locations across methods
                thisVariantInput = [thisVariantInput_border  thisVariantInput_ectopic]; % input 2x locations                
        end
        clear all*vars; % clear to save mem
%%%   NOT WORKING
%     case 'locationXnetwork'
%         load([input_dir 'locationXnetworkInfo.mat']); % starts out as an 814x11 cell, want it to be an 59412 x 11 x 814
%         thisVariantInput = nan*ones(59412,11,814); %initialize
%         switch variant_type
%             case 'border'
%                 for s = 1:size(border_allVars_variantCorrs,1)
%                     for n = 1:size(border_allVars_variantCorrs,2)
%                         thisVariantInput(:,n,s) = border_allVars_variantCorrs{s,n};
%                     end
%                 end
%             case 'ectopic'
%                 for s = 1:size(ectopic_allVars_variantCorrs,1)
%                     for n = 1:size(ectopic_allVars_variantCorrs,2)
%                         thisVariantInput(:,n,s) = ectopic_allVars_variantCorrs{s,n};
%                     end
%                 end
%         end
%         clear *_allVars_variantCorrs; %clear to save mem
%%%
    otherwise
        error('Do not recognize variant form input');
end