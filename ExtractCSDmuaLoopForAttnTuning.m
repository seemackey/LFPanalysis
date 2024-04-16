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

'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004012@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004013@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004014@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004015@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004016@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004017@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004018@os.mat';


};

  paths_vis = {
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004020@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004021@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004022@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004023@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004024@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004025@os.mat';
'C:\Users\cmackey\Documents\AttnTuning\ab004\2-ab003004026@os.mat';


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


    
    [CSD,LFP, MUA, trig, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    outlierfortrigtype = outlieridxtmp <= numel(trigtype);
    trigtype(outlieridxtmp(outlierfortrigtype)) = [];
    % clean data
    data{1,loopct} = num2cell(CSD); 
    data{2,loopct} = num2cell(LFP);
    data{3,loopct} = num2cell(MUA);
    data{4,loopct} = num2cell(trig); % trig times
    data{5,loopct} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    data{6,loopct} = num2cell(std_freq); % tone freq
%     data{1,loopct}(:,outlieridxtmp,:) = []; 
%     data{2,loopct}(:,outlieridxtmp,:) = []; 
%     data{3,loopct}(:,outlieridxtmp,:) = []; 
%     data{4,loopct}(:,outlieridxtmp,:) = []; 
%     data{5,loopct}(outlieridxtmp,:) = [];

    

    
end
toc

%% now we have the data in "data" 
% %% Define channel, time window, and trigger type of interest
selectedChannel = 16:17;  % Change this to the desired channels, ERRORS WHEN USING SINGLE CHANNEL! NOOOOOOO
selectedTimeWindow = 60:85;  % Change this to the desired time window
maximum_previous = 0; %first time we run it we want to calc max, second time we are normalizing attn to ignore, 
%and so we want to norm the second dataset to the max of the first
figure;
subplotArgs = [2 1 1];
plotcolor = 'k';
[maximum,gaussianFit_aud,freqs_aud,avgResp_aud] = AttnTuningFitFxn(data,selectedChannel,selectedTimeWindow,maximum_previous,subplotArgs,plotcolor)

%% extract second data set 

datavis=cell(length(paths_vis));
for loopct = 1:1:length(paths_vis)
    
    [CSD, LFP, trig, trigtype, MUA,epoch_tframe,std_freq] =  EphysExtractFxn(paths_vis{loopct,1},trigch); % extract data and stim triggers


    
    [CSD,LFP, MUA, trig, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    outlierfortrigtype = outlieridxtmp <= numel(trigtype);
    trigtype(outlieridxtmp(outlierfortrigtype)) = [];
    % clean data
    datavis{1,loopct} = num2cell(CSD); 
    datavis{2,loopct} = num2cell(LFP);
    datavis{3,loopct} = num2cell(MUA);
    datavis{4,loopct} = num2cell(trig); % trig times
    datavis{5,loopct} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    datavis{6,loopct} = num2cell(std_freq); % tone freq
%     data{1,loopct}(:,outlieridxtmp,:) = []; 
%     data{2,loopct}(:,outlieridxtmp,:) = []; 
%     data{3,loopct}(:,outlieridxtmp,:) = []; 
%     data{4,loopct}(:,outlieridxtmp,:) = []; 
%     data{5,loopct}(outlieridxtmp,:) = [];


    
end

%% now we have the data 
% %% Defined channel already, time window defined already, and trigger type of interest
maximum_previous = maximum; %first time we run it we want to calc max, 
% second time we are normalizing attn to ignore, and so we want to norm 
% the second dataset to the max of the first
subplotArgs = [2 1 1];
plotcolor = 'r';
[~,gaussianFit_vis,freqs_vis,AvgResps_vis] = AttnTuningFitFxn(datavis,selectedChannel,selectedTimeWindow,maximum_previous,subplotArgs,plotcolor)


% get out parameters and compute difference
% gaussianFitparams_aud = [gaussianFit_aud.a1,gaussianFit_aud.b1,gaussianFit_aud.c1/sqrt(2),gaussianFit_aud.a1+gaussianFit_aud.c1/sqrt(2)];
% 
% gaussianFitparams_vis = [gaussianFit_vis.a1,gaussianFit_vis.b1,gaussianFit_vis.c1/sqrt(2),gaussianFit_vis.a1+gaussianFit_vis.c1/sqrt(2)];
% 
% ratios = gaussianFitparams_aud./gaussianFitparams_vis;

coeffs_aud = coeffvalues(gaussianFit_aud);
gaussianFitparams_aud = [coeffs_aud(1),coeffs_aud(2),coeffs_aud(3)/sqrt(2),coeffs_aud(1)+coeffs_aud(3)/sqrt(2)];

coeffs_vis = coeffvalues(gaussianFit_vis);
gaussianFitparams_vis = [coeffs_vis(1),coeffs_vis(2),coeffs_vis(3)/sqrt(2),coeffs_vis(1)+coeffs_vis(3)/sqrt(2)];

ratios = gaussianFitparams_aud./gaussianFitparams_vis;

% %%
%clear data datavis plotcolor trig LFP CSD MUA maximum maximum_previous gimmeplots 

%% Save the workspace
%workspace_filename = strcat(paths{1,1}, '_workspace.mat');
%save(workspace_filename,'-v7.3');

% Save the current figure
%fig_filename = strcat(paths{1,1}, '_figure.fig');  % Modify as needed
%saveas(gcf, fig_filename);