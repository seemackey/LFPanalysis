% This loops through the ephys extract fxn and gets MUA out, then checks 
% if the response is significantly greater than baseline

clear
close all

%% set directory
myDir = 'G:\A1_LED files'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct


%% loop through directory
for loopct = 1:length(myFiles)
    
    baseFileName = myFiles(loopct).name;
    fullfilename = fullfile(myDir, baseFileName);
  
    [CSD, LFP, trig, trigtype, MUA] =  EphysExtractFxn(fullfilename,3); % extract data and stim triggers

    %unclean data
    all{1,1} = num2cell(CSD); 
    all{2,1} = num2cell(LFP);
    all{3,1} = num2cell(MUA);
    all{4,1} = num2cell(trig); % trig times
    all{5,1} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    
    [~, ~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    
    
    % clean data
    all{1,1}(:,outlieridxtmp,:) = []; 
    all{2,1}(:,outlieridxtmp,:) = []; 
    all{3,1}(:,outlieridxtmp,:) = []; 
    all{4,1}(:,outlieridxtmp,:) = []; 
    all{5,1}(outlieridxtmp,:) = [];
    
    tpmua1=round(length(MUA(1,1,:))/2); %using halfway point because I've been epoching 250 ms before stim and 250 after
    tpmua2=length(MUA(1,1,:));

    MUA = cell2mat(all{3,1}(:,:,:)); % MUA 
    trigtype = cell2mat(all{5,1}(:,:)); 
    
    for i=1:size(MUA,1)
        for iii=1:size(MUA,2)
            for times_smooth=1:2
             MUA(i,iii,:)=smooth(squeeze(MUA(i,iii,:)),'moving',5);%smoothed data
            end
        end
    end
    
%% calculate confidence intervals of MUA and significance of response
    numchans = length(all{1,1}(:,1,1));
    numtrs = length(all{1,1}(1,:,1));
    numtps = length(squeeze(all{1,1}(1,1,:)));
    ci_stp = 1; % step size, multiply by 2 to get milliseconds
    numboots = 1000; % number of resamples
    
    sigtimes_sum = zeros(length(numchans),tpmua2);
    % loop through channels
    for ci_ch_ct = 1:1:numchans
        for ci_tp = 1:ci_stp:tpmua2 % loop through time points
            MUAtmp = squeeze(MUA(ci_ch_ct,:,ci_tp)); % get the data from a time point
            ci_tmp = bootci(numboots,@mean,MUAtmp'); % get bootstrap derived CIs
            cimean(ci_ch_ct,ci_tp) = mean(MUAtmp); % record mean
            cilo(ci_ch_ct,ci_tp) = ci_tmp(1,1); % lower CI
            cihi(ci_ch_ct,ci_tp) = ci_tmp(2,1); % upper CI

            if cilo(ci_ch_ct,ci_tp)>0
                sigtimes(ci_ch_ct,ci_tp) = 1;
            elseif cihi(ci_ch_ct,ci_tp)<0
                sigtimes(ci_ch_ct,ci_tp) = -1;
            else
                sigtimes(ci_ch_ct,ci_tp) = 0;
            end
        end
    end
%% % add sig times
  sigtimes_all{loopct} = sigtimes;
  sigtimes_sum = sigtimes+sigtimes_sum;
  
  
  cihi_all{loopct,1} = cihi;
  cilo_all{loopct,1} = cilo;
  cimeanall{loopct,1} = cimean;
  
  
end
save('sigtimes');