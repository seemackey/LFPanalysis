%% Initialization
tic;
clear;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\AttnData\core\cont\Figures\'; % Specify directory
myfiles = dir(fullfile(mydir, '*AveragePowerPhaseData.mat'));  % Adjust the pattern to match your new data files

%% Process each file
% Initialize arrays to store the channels with the highest average power
highestChannelCSD = [];
highestChannelLFP = [];

for loopct = 1:1%length(myfiles)
    % Load the data file
    dataFilePath = fullfile(mydir, myfiles(loopct).name);
    load(dataFilePath,'avgPhaseCSDDesc','avgPhaseLFPDesc', 'avgPowerCSDDesc', 'avgPowerLFPDesc', 'channelsDesc');
    
    % Find the indices of the maximum average power for CSD and LFP
    [~, idxCSD] = max(avgPowerCSDDesc);
    [~, idxLFP] = max(avgPowerLFPDesc);
    
    % Store the channel with the highest average power
    highestChannelCSD = [highestChannelCSD; channelsDesc(idxCSD)];
    highestChannelLFP = [highestChannelLFP; channelsDesc(idxLFP)];
end

%% Create histograms (channel counts)
channelCountCSD = histcounts(highestChannelCSD, 0.5:1:(max(channelsDesc)+0.5));
channelCountLFP = histcounts(highestChannelLFP, 0.5:1:(max(channelsDesc)+0.5));

%% Plotting histograms as horizontal bar graphs
figure;
% Horizontal bar graph for CSD
subplot(1, 2, 1);
barh(1:length(channelCountCSD), channelCountCSD, 'FaceColor', 'b');
title('Max. Average Power (CSD)');
xlabel('Count');
ylabel('Channel');
set(gca, 'YDir', 'reverse', 'YTick', 1:length(channelsDesc), 'YTickLabel', flip(channelsDesc));  % Ensure correct channel labeling in descending order
axis tight;

% Horizontal bar graph for LFP
subplot(1, 2, 2);
barh(1:length(channelCountLFP), channelCountLFP, 'FaceColor', 'r');
title('Max. Average Power (LFP)');
xlabel('Count');
ylabel('Channel');
set(gca, 'YDir', 'reverse', 'YTick', 1:length(channelsDesc), 'YTickLabel', flip(channelsDesc));  % Ensure correct channel labeling in descending order
axis tight;

%% Finalization
toc;
