%% Loops through the ephys extract function and gets out LFP data

clear
close all

% make this 1 if you want to see acoustic responses 
gimmeplots=1;

% select trigger, 1 for aud, 3 for visual, 4 for sacc on, 5 for sacc off
trigch=1;

%% path info 

%  paths = {'D:\oldtono\contproc\1-gt030031021@os_eye06_30'
% 'D:\oldtono\contproc\1-gt030031021@os_eye06_30'};

%  paths = {
%      '/Volumes/14TB_USB_01/jitter/_missing/ma041/2-ma041042041@os.mat';
% };

  paths = {
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042017@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042018@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042019@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042020@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042021@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042022@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042023@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042024@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042025@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042026@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042027@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042028@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042029@os';
};

  paths_vis = {
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042033@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042034@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042035@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042036@os';
'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042037@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042038@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042039@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042040@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042041@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042042@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042043@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042044@os';
 'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042045@os';
  'C:\Users\cmackey\Documents\AttnTuning\kk041\1-kk041042046@os';
};
 
 %/Volumes/16TB_003/dyneyep/bbn/contproc/2-rb069070030@os_eye06_20.mat bbn

% mac path example
% paths = {'/Volumes/Samsung03/bbn/1-bu015016038@os.mat'
%     '/Volumes/Samsung03/bbn/2-bu015016038@os.mat'};


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\1-bu019020039@os_eye06_20.mat';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\2-bu019020039@os_eye06_20.mat'; };

%% loop through the function to extract data and find outliers
data=cell(length(paths));
tic
for loopct = 1:1:length(paths)
    
    [CSD, LFP, trig, trigtype, MUA,epoch_tframe,std_freq] =  EphysExtractFxn(paths{loopct,1},trigch); % extract data and stim triggers

    %unclean data
    data{1,loopct} = num2cell(CSD); 
    data{2,loopct} = num2cell(LFP);
    data{3,loopct} = num2cell(MUA);
    data{4,loopct} = num2cell(trig); % trig times
    data{5,loopct} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    data{6,loopct} = num2cell(std_freq); % tone freq
    
    [~, ~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    
    data{1,outlierct}(:,outlieridxtmp,:) = []; 
    data{2,outlierct}(:,outlieridxtmp,:) = []; 
    data{3,outlierct}(:,outlieridxtmp,:) = []; 
    data{4,outlierct}(:,outlieridxtmp,:) = []; 
    data{5,outlierct}(outlieridxtmp,:) = [];
    data{6,outlierct}(outlieridxtmp,:) = [];
    
    outliers{:,loopct} = outlieridxtmp; %at the end this has idx of outliers from both recordings

    
end
toc

%% now we have the data in "data" 
% %% Define channel, time window, and trigger type of interest
selectedChannel = 12:16;  % Change this to the desired channel
selectedTimeWindow = 55:90;  % Change this to the desired time window
maximum_previous = 0; %first time we run it we want to calc max, second time we are normalizing attn to ignore, 
%and so we want to norm the second dataset to the max of the first
figure;
subplotArgs = [2 1 1];
plotcolor = 'k';
[maximum,gaussianFit_aud,freqs_aud] = AttnTuningFitFxn(data,selectedChannel,selectedTimeWindow,maximum_previous,subplotArgs,plotcolor)

%% extract second data set 

datavis=cell(length(paths_vis));
for loopct = 1:1:length(paths_vis)
    
    [CSD, LFP, trig, trigtype, MUA,epoch_tframe,std_freq] =  EphysExtractFxn(paths_vis{loopct,1},trigch); % extract data and stim triggers

    %unclean data
    datavis{1,loopct} = num2cell(CSD); 
    datavis{2,loopct} = num2cell(LFP);
    datavis{3,loopct} = num2cell(MUA);
    datavis{4,loopct} = num2cell(trig); % trig times
    datavis{5,loopct} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    datavis{6,loopct} = num2cell(std_freq); % tone freq
    
    [~, ~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    
    datavis{1,outlierct}(:,outlieridxtmp,:) = []; 
    datavis{2,outlierct}(:,outlieridxtmp,:) = []; 
    datavis{3,outlierct}(:,outlieridxtmp,:) = []; 
    datavis{4,outlierct}(:,outlieridxtmp,:) = []; 
    datavis{5,outlierct}(outlieridxtmp,:) = [];
    datavis{6,outlierct}(outlieridxtmp,:) = [];
    [~, ~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    
    outliers{:,loopct} = outlieridxtmp; %at the end this has idx of outliers from both recordings
end

%% now we have the data in "all"
% %% Defined channel already, time window defined already, and trigger type of interest
maximum_previous = maximum; %first time we run it we want to calc max, 
% second time we are normalizing attn to ignore, and so we want to norm 
% the second dataset to the max of the first
subplotArgs = [2 1 1];
plotcolor = 'r';
[~,gaussianFit_vis,freqs_vis] = AttnTuningFitFxn(datavis,selectedChannel,selectedTimeWindow,maximum_previous,subplotArgs,plotcolor)


%% get out parameters and compute difference
gaussianFitparams_aud = [gaussianFit_aud.a1,gaussianFit_aud.b1,gaussianFit_aud.c1/sqrt(2),gaussianFit_aud.a1+gaussianFit_aud.c1/sqrt(2)];

gaussianFitparams_vis = [gaussianFit_vis.a1,gaussianFit_vis.b1,gaussianFit_vis.c1/sqrt(2),gaussianFit_vis.a1+gaussianFit_vis.c1/sqrt(2)];

ratios = gaussianFitparams_aud./gaussianFitparams_vis;

clear plotcolor trig LFP CSD MUA maximum maximum_previous gimmeplots 