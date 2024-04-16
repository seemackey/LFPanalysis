% find peak power in wavelet files

%% Initialization
tic
clear
close all

%% Set directory and parameters
mydir = 'H:\jitter\gt038_02\'; % Specify directory
myfiles = dir(fullfile(mydir,'*@osw*')); % Get all files in struct

freq_range = [8, 14]; % Frequency range of interest (Hz)
pre_stim_ms = -300; % Time window start relative to stimulus (ms)
post_stim_ms = -100; % Time window end relative to stimulus (ms)
%samplingRate = 1000; % Example sampling rate in Hz, replace with wraw.adrate as necessary

%% Data collection for heatmap
allPeakTimes = []; % Collect all peak times across files
allChannelIdx = []; % Collect corresponding channel indices for peaks

%% Process each file
for loopct = 1:1%length(myfiles)
    
    %% Load the wavelet file
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    load(fullfilename); % 
    samplingRate = wraw.adrate;
    %% Identify trigger times of type "1"
    triggerTypeOnesIndices = find(trig.ttype{1,1} == 1);
    triggerTimesTypeOne = trig.anatrig{1,1}(triggerTypeOnesIndices);

    %% Process each trigger
    for i = 1:length(triggerTimesTypeOne)
        triggerTime = triggerTimesTypeOne(i);

        % Define epoch window around the trigger
        epochStart = triggerTime + round(pre_stim_ms * (samplingRate / 1000));
        epochEnd = triggerTime + round(post_stim_ms * (samplingRate / 1000));
        
        % Ensure epoch is within data bounds
        epochStart = max(epochStart, 1);
        epochEnd = min(epochEnd, size(wraw.cnte_po, 3));

        % Extract epoch data for the frequency range at curr. trig.
        freq_indices = find(wraw.frq >= freq_range(1) & wraw.frq <= freq_range(2));
        epochDataCSD = wraw.cntc_po(:, freq_indices, epochStart:epochEnd);
        epochDataLFP = wraw.cnte_po(:, freq_indices, epochStart:epochEnd);
        epochDataCSDph = wraw.cntc_ph(:, freq_indices, epochStart:epochEnd);
        epochDataLFPph = wraw.cnte_ph(:, freq_indices, epochStart:epochEnd);

        % CSD    
        % Find peak power within epoch
        [peakPowerCSD, idxCSD] = max(epochDataCSD(:));
        [chIdxCSD, ~, tIdxCSD] = ind2sub(size(epochDataCSD), idxCSD);
        
        % Adjust peak time relative to trigger and convert to ms
        peakTime = ((epochStart + tIdxCSD - 1) - triggerTime) * (1000 / samplingRate);
        
        %LFP
        % Find peak power within epoch
        [peakPowerLFP, idxLFP] = max(epochDataLFP(:));
        [chIdxLFP, ~, tIdxLFP] = ind2sub(size(epochDataLFP), idxLFP);
        
        % Adjust peak time relative to trigger and convert to ms
        peakTime = ((epochStart + tIdxCSD - 1) - triggerTime) * (1000 / samplingRate);
        
        %phases
        % Get phases of power peaks using ind2sub
        [chIdxCSDph,~,tIdxCSDph] = ind2sub(size(epochDataCSDph),idxCSD); % does this make sense?
       
        
        epochDataLFPph = []; % same for FP
        
        % Collect data for heatmap
        allPowerPeaksCSD = [allPeakPowerPeaksCSD;peakPowerCSD];
        allPeakTimesCSD = [allPeakTimesCSD; peakTimeCSD];
        allChannelIdxCSD = [allChannelIdxCSD; chIdxCSD];
    end
end

%% Generate heatmap of peak times
% Define bins for heatmap
% Adjust timeEdges for 2 ms bins from -300 ms to 100 ms
timeEdges = linspace(-300, 10, 21); % Creates 200 bins of 2 ms each

% Define channelEdges for one bin per channel, ensuring channels are plotted descending
channelEdges = 0.5:1:(max(allChannelIdx)+0.5); % One bin per channel

% Create 2D histogram of peak times vs. channels
[peakCounts, ~, ~] = histcounts2(allChannelIdx, allPeakTimes, channelEdges, timeEdges);
% Flip the peakCounts matrix vertically so channel 1's data appears at the top of the heatmap
%peakCountsFlipped = flipud(peakCounts);

% Plot heatmap
figure;
imagesc('XData', timeEdges(1:end-1), 'YData', 1:max(allChannelIdx), 'CData', peakCounts);
colormap('jet'); % Set the colormap to 'jet'
colorbar; % Show color scale
xlabel('Time Relative to Trigger (ms)');
ylabel('Channel');
title('Distribution of Peak 8-14 Hz Activity Across Channels and Time');
axis tight; % Fit the axes to the data
set(gca, 'YDir', 'reverse');




%% Finalization

% Assuming processing for one of the files is done and you're at the saving step

% Base file name for saving
processedFileName = basefilename(1:end-4); % Remove file extension if present

% Define file paths for saving
variablesFileName = fullfile(mydir, [processedFileName, '_ProcessedVariables.mat']);
figureFileName = fullfile(mydir, [processedFileName, '_PeakActivityHeatmap.png']);

% Save key variables
save(variablesFileName, 'allPeakTimes', 'allChannelIdx', 'peakCounts', 'timeEdges', 'channelEdges');

% Save figure
figureHandle = figure;
imagesc('XData', timeEdges(1:end-1), 'YData', 1:max(allChannelIdx), 'CData', peakCounts);
colormap('jet');
colorbar;
xlabel('Time Relative to Trigger (ms)');
ylabel('Channel');
title('Distribution of Peak 8-14 Hz Activity Across Channels and Time');
axis tight;
xticks(linspace(-300, 300, 13));
set(gca, 'YDir', 'reverse');

% For MATLAB versions before R2019b, use saveas:
saveas(figureHandle, figureFileName);

% Close the figure if you don't want it to stay open
close(figureHandle);

toc
