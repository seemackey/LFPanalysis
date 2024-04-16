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

% Plot the first dataset
subplotArgs1 = [2 1 1];
plotcolor1 = 'k';
[maximum1, gaussianFit_aud, freqs_aud] = AttnTuningFitFxnv2(data, selectedChannel, selectedTimeWindow, 0, subplotArgs1, plotcolor1);

% Extract second dataset
clear datavis
datavis = cell(length(paths_vis), 6);
for loopct = 1:length(paths_vis)
    [CSD, LFP, trig, trigtype, MUA, epoch_tframe, std_freq] = EphysExtractFxn(paths_vis{loopct, 1}, trigch);
    datavis{loopct, 1} = CSD;
    datavis{loopct, 2} = LFP;
    datavis{loopct, 3} = MUA;
    datavis{loopct, 4} = trig;
    datavis{loopct, 5} = trigtype;
    datavis{loopct, 6} = std_freq;
end

% Plot the second dataset
subplotArgs2 = [2 1 2];
plotcolor2 = 'r';
[maximum2, gaussianFit_vis, freqs_vis] = AttnTuningFitFxnv2(datavis, selectedChannel, selectedTimeWindow, maximum1, subplotArgs2, plotcolor2);

% Clear the second dataset variable
clear datavis;