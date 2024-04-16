% snr wavelet avg power
% loops through wavelet files and gets avg (over trials and time) power at
% the freq range specified below
%% Initialization
tic;
clear;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\AttnData\core\wavelet\'; % Specify directory
myfiles = dir(fullfile(mydir,'*@osw*')); % Get all files in struct

freq_range = [8, 14]; % Frequency range of interest (Hz)
pre_stim_ms = -300; % Time window start relative to stimulus (ms)
post_stim_ms = -100; % Time window end relative to stimulus (ms)

%% Initialize variables to store power values across all files and triggers
avgPowerCSD = [];
avgPowerLFP = [];
stdPowerCSD = [];
stdPowerLFP = [];
numChannels = 0;

%% Process each file
for loopct = 1:length(myfiles)
    
    %% Load the wavelet file
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    load(fullfilename); % Load your .mat file containing 'wraw', 'trig', etc.
    samplingRate = wraw.adrate; % Assuming this is defined within each loaded file
    
    %% Identify trigger times of type "1"
    triggerTypeOnesIndices = find(trig.ttype{1,1} == 1);
    triggerTimesTypeOne = trig.anatrig{1,1}(triggerTypeOnesIndices);

    %% Initialize temporary storage for this file's calculations
    tempPowerCSD = zeros(length(triggerTimesTypeOne), size(wraw.cntc_po, 1));
    tempPowerLFP = zeros(length(triggerTimesTypeOne), size(wraw.cnte_po, 1));
    
    %% Process each trigger
    tempAvgPowerCSD = zeros(length(triggerTimesTypeOne), size(wraw.cntc_po, 1));
    tempAvgPowerLFP = zeros(length(triggerTimesTypeOne), size(wraw.cnte_po, 1));

    for i = 1:length(triggerTimesTypeOne)
        triggerTime = triggerTimesTypeOne(i);

        % Define epoch window around the trigger
        epochStart = triggerTime + round(pre_stim_ms * (samplingRate / 1000));
        epochEnd = triggerTime + round(post_stim_ms * (samplingRate / 1000));

        % Ensure epoch is within data bounds
        epochStart = max(epochStart, 1);
        epochEnd = min(epochEnd, size(wraw.cnte_po, 3));

        % Extract epoch data for the frequency range
        freq_indices = find(wraw.frq >= freq_range(1) & wraw.frq <= freq_range(2));
        epochDataCSD = wraw.cntc_po(:, freq_indices, epochStart:epochEnd);
        epochDataLFP = wraw.cnte_po(:, freq_indices, epochStart:epochEnd);

        % Average power over time and frequency for the epoch
        % Each row in tempAvgPowerCSD and tempAvgPowerLFP corresponds to a trial,
        % and each column corresponds to a channel.
        tempAvgPowerCSD(i, :) = mean(mean(epochDataCSD, 3), 2)'; % Average over time, then frequency
        tempAvgPowerLFP(i, :) = mean(mean(epochDataLFP, 3), 2)'; % Average over time, then frequency
    end

    % After processing all triggers, compute the average and standard deviation across trials (rows of tempAvgPower*)
    avgPowerCSD = mean(tempAvgPowerCSD, 1); % Average across trials for each channel
    avgPowerLFP = mean(tempAvgPowerLFP, 1); % Average across trials for each channel
    stdPowerCSD = std(tempAvgPowerCSD, 0, 1); % STD across trials for each channel
    stdPowerLFP = std(tempAvgPowerLFP, 0, 1); % STD across trials for each channel

    
    % Update number of channels if necessary
    numChannels = size(wraw.cntc_po, 1);
end

%% Plotting the average power for each electrode channel with error bars

% Prepare the channel numbers as y-axis values
channels = 2:numChannels; % Assuming numChannels was set to the number of channels earlier
channelsDesc = flip(channels); % Flip the channel numbers for descending order plotting

% Preparing data for descending channel order
avgPowerCSDDesc = avgPowerCSD(channelsDesc);
avgPowerLFPDesc = avgPowerLFP(channelsDesc);
stdPowerCSDDesc = stdPowerCSD(channelsDesc);
stdPowerLFPDesc = stdPowerLFP(channelsDesc);

figure;

% CSD Power
subplot(1, 2, 1);
errorbar(avgPowerCSDDesc, channels, stdPowerCSDDesc, 'horizontal');
title('Average CSD Power (8-14 Hz)');
xlabel('Power');
ylabel('Channel');
%set(gca, 'YDir', 'reverse'); % Reverse the y-axis to have channel 1 at the top
set(gca, 'YTick', channels, 'YTickLabel', channelsDesc); % Ensure correct channel labeling
ylim([1 numChannels]);

% LFP Power
subplot(1, 2, 2);
errorbar(avgPowerLFPDesc, channels, stdPowerLFPDesc, 'horizontal');
title('Average LFP Power (8-14 Hz)');
xlabel('Power');
ylabel('Channel');
%set(gca, 'YDir', 'reverse'); % Reverse the y-axis to have channel 1 at the top
set(gca, 'YTick', channels, 'YTickLabel', channelsDesc); % Ensure correct channel labeling
ylim([1 numChannels]);


%% After computing avgPowerCSDDesc, avgPowerLFPDesc, stdPowerCSDDesc, stdPowerLFPDesc

% Define the filenames for saving the figures and data
figureFileNameCSD = fullfile(mydir, [basefilename(1:end-4), '_AveragePower_CSD.fig']);
figureFileNameLFP = fullfile(mydir, [basefilename(1:end-4), '_AveragePower_LFP.fig']);
dataFileName = fullfile(mydir, [basefilename(1:end-4), '_AveragePowerData.mat']);

% Plot and save CSD Power figure
fCSD = figure;
errorbar(avgPowerCSDDesc, channels, stdPowerCSDDesc, 'horizontal');
title('Average CSD Power (8-14 Hz)');
xlabel('Power');
ylabel('Channel');
%set(gca, 'YDir', 'reverse');
set(gca, 'YTick', channels, 'YTickLabel', channelsDesc);
ylim([1 numChannels]);
saveas(fCSD, figureFileNameCSD);

% Plot and save LFP Power figure
fLFP = figure;
errorbar(avgPowerLFPDesc, channels, stdPowerLFPDesc, 'horizontal');
title('Average LFP Power (8-14 Hz)');
xlabel('Power');
ylabel('Channel');
%set(gca, 'YDir', 'reverse');
set(gca, 'YTick', channels, 'YTickLabel', channelsDesc);
ylim([1 numChannels]);
saveas(fLFP, figureFileNameLFP);

% Save the data
save(dataFileName, 'avgPowerCSDDesc', 'avgPowerLFPDesc', 'stdPowerCSDDesc', 'stdPowerLFPDesc', 'channelsDesc');
