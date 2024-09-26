%% Initialization
clc;
clear;
close all;

% Directory settings
imgdir = 'E:\spectrolaminar\AttnData\core\cont_offBF\';
datadir = 'E:\spectrolaminar\AttnData\core\cont\GammaPost50\';

% Get list of image and data files
imgFiles = dir(fullfile(imgdir, '*.jpg'));  % Assuming image files are in jpg format
dataFiles = dir(fullfile(datadir, '*@os_AveragePowerPhaseData.mat'));

% Preallocate list for selected channels
selectedChannels = zeros(length(dataFiles), 1);

%% Loop to display images and collect channel numbers
for k = 1:length(imgFiles)
    % Display image
    figure;
    img = imread(fullfile(imgdir, imgFiles(k).name));
    imshow(img);
    title(['Select the channel for file: ', imgFiles(k).name]);
    
    % User input for channel
    selectedChannels(k) = input('Enter the channel number aligned with the feature: ');
    close all;
end

%%
% Preallocate big array for averaging, size depending on the number of channels found in the first file
% Determine the reference channel for alignment based on the first file
referenceChannel = selectedChannels(1);

% Assess the maximum number of channels in any dataset to set array bounds
maxChannels = 0;
for k = 1:length(dataFiles)
    data = load(fullfile(datadir, dataFiles(k).name), 'avgPowerCSDDesc');
    maxChannels = max(maxChannels, length(data.avgPowerCSDDesc));
end

% Set total rows in the aligned arrays to accommodate the maximum shift
totalRows = 2 * maxChannels;  % Enough space to avoid data loss

% Initialize the big arrays with NaN
alignedPowerCSD = NaN(totalRows, length(dataFiles));
alignedPowerLFP = NaN(totalRows, length(dataFiles));

%% Align and store data for each file
for k = 1:length(dataFiles)
    data = load(fullfile(datadir, dataFiles(k).name), 'avgPowerCSDDesc', 'avgPowerLFPDesc');
    currentChannels = length(data.avgPowerCSDDesc);

    % Calculate the shift needed to align the selected channel with the reference
    shiftIndex = selectedChannels(k) - referenceChannel;
    centerIdx = totalRows / 2;  % Middle of the big array

    % Calculate start and end indices for the data insertion
    startIdx = centerIdx - shiftIndex;
    endIdx = startIdx + currentChannels - 1;

    % Insert the data into the big array using calculated indices
    alignedPowerCSD(startIdx:endIdx, k) = data.avgPowerCSDDesc;
    alignedPowerLFP(startIdx:endIdx, k) = data.avgPowerLFPDesc;
end

%% Calculate averages and SEM ignoring NaNs

alignedPowerCSD(alignedPowerCSD == 0) = NaN;
alignedPowerLFP(alignedPowerLFP == 0) = NaN;


avgPowerCSD = nanmean(alignedPowerCSD, 2);
semPowerCSD = nanstd(alignedPowerCSD, 0, 2) / sqrt(sum(~isnan(alignedPowerCSD), 2));
semPowerCSD(semPowerCSD == 0) = NaN;

avgPowerLFP = nanmean(alignedPowerLFP, 2);
semPowerLFP = nanstd(alignedPowerLFP, 0, 2) / sqrt(sum(~isnan(alignedPowerLFP), 2));
semPowerLFP(semPowerLFP == 0) = NaN;

figure;
subplot(1,2,1);
hold on;
errorbar(alignedXCSD, avgPowerCSD(validCSD), semPowerCSD(~isnan(semPowerCSD)), 'horizontal', 'k', 'LineWidth', 2);
title('CSD Gamma Amp.');
xlabel('Distance from Inversion (mm)');
ylabel('Normalized Power');
set(gca, 'YDir', 'reverse');  % Ensure the highest channel number is at the top of the plot
hold off;

subplot(1,2,2);
hold on;
errorbar(alignedXLFP, avgPowerLFP(validLFP), semPowerLFP(~isnan(semPowerLFP)), 'horizontal', 'k', 'LineWidth', 2);
title('LFP Gamma Amp.');
xlabel('Distance from Inversion (mm)');
ylabel('Normalized Power');
set(gca, 'YDir', 'reverse');  % Ensure the highest channel number is at the top of the plot
hold off;

%%
% Parameters for bootstrapping
numBootstraps = 1000;  % Number of bootstrap samples
bootstrapFunc = @(x) nanmean(x, 2);  % Function to compute mean, ignoring NaNs

% Initialize the CI arrays
ciLowerCSD = NaN(size(alignedPowerCSD, 1), 1);
ciUpperCSD = NaN(size(alignedPowerCSD, 1), 1);
ciLowerLFP = NaN(size(alignedPowerLFP, 1), 1);
ciUpperLFP = NaN(size(alignedPowerLFP, 1), 1);

% Minimum number of non-NaN observations required
minNonNaN = 3;

% Calculate CIs for each channel, checking for sufficient data
for i = 1:size(alignedPowerCSD, 1)
    if sum(~isnan(alignedPowerCSD(i, :))) >= minNonNaN
        ci = bootci(numBootstraps, {bootstrapFunc, alignedPowerCSD(i, :)}, 'type', 'per');
        ciLowerCSD(i) = ci(1);
        ciUpperCSD(i) = ci(2);
    end

    if sum(~isnan(alignedPowerLFP(i, :))) >= minNonNaN
        ci = bootci(numBootstraps, {bootstrapFunc, alignedPowerLFP(i, :)}, 'type', 'per');
        ciLowerLFP(i) = ci(1);
        ciUpperLFP(i) = ci(2);
    end
end
%%
% Assuming that avgPowerCSD and avgPowerLFP have been computed
% Calculate the indices of non-NaN entries for both CSD and LFP
validCSD = ~isnan(avgPowerCSD) & sum(~isnan(alignedPowerCSD), 2) >= 3;
validLFP = ~isnan(avgPowerLFP) & sum(~isnan(alignedPowerLFP), 2) >= 3;

% Calculate x-axis values only for valid data
alignedXCSD = 0.1 * ((1:sum(validCSD)) - selectedChannels(1));  % Adjusted for valid CSD indices
alignedXLFP = 0.1 * ((1:sum(validLFP)) - selectedChannels(1));  % Adjusted for valid LFP indices

figure;
subplot(1,2,1);
hold on;
% Shading for CSD
fill([alignedXCSD, fliplr(alignedXCSD)], [ciLowerCSD(validCSD); flipud(ciUpperCSD(validCSD))], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(alignedXCSD, avgPowerCSD(validCSD), 'k', 'LineWidth', 2); % Avg power CSD
title('Aligned CSD Amp.');
xlabel('Distance from inversion (mm)');
ylabel('Average Power');
hold off;

subplot(1,2,2);
hold on;
% Shading for LFP
fill([alignedXLFP, fliplr(alignedXLFP)], [ciLowerLFP(validLFP); flipud(ciUpperLFP(validLFP))], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(alignedXLFP, avgPowerLFP(validLFP), 'k', 'LineWidth', 2); % Avg power LFP
title('Aligned LFP Amp.');
xlabel('Distance from inversion (mm)');
ylabel('Average Power');
hold off;


%% Plotting results for CSD
clean = ~isnan(avgPowerCSD);
alignedX = alignedXCSD;

%
% figure
% subplot(1,2,1)
% cleanSEM = ~isnan(semPowerCSD);
% plot(avgPowerCSD(clean), alignedX')
% title('Aligned CSD Alpha Post Amp.');
% ylabel('Distance from Inversion (mm)');
% subplot(1,2,2)
% cleanSEM = ~isnan(semPowerLFP);
% plot(avgPowerLFP(clean), alignedX')
% title('Aligned LFP Alpha Post Amp.');
% ylabel('Distance from Inversion (mm)');

figure
subplot(1,2,1)
plot(avgPowerCSD(validCSD),alignedXCSD,  'k', 'LineWidth', 2); % Avg power CSD
title('CSD Gamma Amp.');
ylabel('Distance from Inversion (mm)');
subplot(1,2,2)
plot(avgPowerLFP(validLFP),alignedXLFP,  'k', 'LineWidth', 2); % Avg power LFP
title('LFP Gamma Amp.');
ylabel('Distance from Inversion (mm)');

%% NORMALIZE AND BOOTSTRAP
% Initialize the normalized power arrays
normalizedPowerCSD = NaN(size(alignedPowerCSD));
normalizedPowerLFP = NaN(size(alignedPowerLFP));

% Normalize each column to its own maximum channel power
for col = 1:size(alignedPowerCSD, 2)
    % Find the maximum power in each column for CSD
    maxPowerColCSD = max(alignedPowerCSD(:, col));
    if ~isnan(maxPowerColCSD) && maxPowerColCSD ~= 0
        normalizedPowerCSD(:, col) = alignedPowerCSD(:, col) / maxPowerColCSD;
    end

    % Find the maximum power in each column for LFP
    maxPowerColLFP = max(alignedPowerLFP(:, col));
    if ~isnan(maxPowerColLFP) && maxPowerColLFP ~= 0
        normalizedPowerLFP(:, col) = alignedPowerLFP(:, col) / maxPowerColLFP;
    end
end

% Recalculate average power for normalized data
avgNormalizedPowerCSD = nanmean(normalizedPowerCSD, 2);
avgNormalizedPowerLFP = nanmean(normalizedPowerLFP, 2);


% Perform bootstrapping for normalized data
% Parameters for bootstrapping
numBootstraps = 1000;  % Number of bootstrap samples
bootstrapFunc = @(x) nanmean(x, 2);  % Function to compute mean, ignoring NaNs

% Initialize the CI arrays for normalized data
ciLowerNormCSD = NaN(size(normalizedPowerCSD, 1), 1);
ciUpperNormCSD = NaN(size(normalizedPowerCSD, 1), 1);
ciLowerNormLFP = NaN(size(normalizedPowerLFP, 1), 1);
ciUpperNormLFP = NaN(size(normalizedPowerLFP, 1), 1);

% Calculate CIs for each channel, checking for sufficient data
for i = 1:size(normalizedPowerCSD, 1)
    if sum(~isnan(normalizedPowerCSD(i, :))) >= minNonNaN
        ci = bootci(numBootstraps, {bootstrapFunc, normalizedPowerCSD(i, :)}, 'type', 'per');
        ciLowerNormCSD(i) = ci(1);
        ciUpperNormCSD(i) = ci(2);
    end

    if sum(~isnan(normalizedPowerLFP(i, :))) >= minNonNaN
        ci = bootci(numBootstraps, {bootstrapFunc, normalizedPowerLFP(i, :)}, 'type', 'per');
        ciLowerNormLFP(i) = ci(1);
        ciUpperNormLFP(i) = ci(2);
    end
end

% Plotting normalized power with confidence intervals
figure;
subplot(1,2,1);
hold on;
validNormCSD = ~isnan(avgNormalizedPowerCSD) & sum(~isnan(normalizedPowerCSD), 2) >= 3;
fill([alignedXCSD, fliplr(alignedXCSD)], [ciLowerNormCSD(validNormCSD); flipud(ciUpperNormCSD(validNormCSD))], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(alignedXCSD, avgNormalizedPowerCSD(validNormCSD), 'k', 'LineWidth', 2); % Normalized power CSD
title('Normalized Aligned CSD Alpha Amp.');
xlabel('Distance from inversion (mm)');
ylabel('Normalized Average Power');
hold off;

subplot(1,2,2);
hold on;
validNormLFP = ~isnan(avgNormalizedPowerLFP) & sum(~isnan(normalizedPowerLFP), 2) >= 5;
fill([alignedXLFP, fliplr(alignedXLFP)], [ciLowerNormLFP(validNormLFP); flipud(ciUpperNormLFP(validNormLFP))], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(alignedXLFP, avgNormalizedPowerLFP(validNormLFP), 'k', 'LineWidth', 2); % Normalized power LFP
title('Normalized Aligned LFP Alpha Amp.');
xlabel('Distance from inversion (mm)');
ylabel('Normalized Average Power');
hold off;

%%
figure;

% % Calculate valid indices for CSD and LFP
% validNormCSD = ~isnan(avgNormalizedPowerCSD) & sum(~isnan(normalizedPowerCSD), 2) >= 5;
% validNormLFP = ~isnan(avgNormalizedPowerLFP) & sum(~isnan(normalizedPowerLFP), 2) >= 5;
% 
% % Define x-axis for the channels based on the valid indices calculated above
% alignedXCSD = 0.1 * ((1:length(avgNormalizedPowerCSD)) - selectedChannels(1));
% alignedXLFP = 0.1 * ((1:length(avgNormalizedPowerLFP)) - selectedChannels(1));
% 
% % Filter x-axis values based on valid normalized indices
% alignedXCSD = alignedXCSD(validNormCSD);
% alignedXLFP = alignedXLFP(validNormLFP);

% Plotting for CSD
subplot(1,2,1);
hold on;
fill([ciLowerNormCSD(validNormCSD); flipud(ciUpperNormCSD(validNormCSD))], ...
     [alignedXCSD, fliplr(alignedXCSD)], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(avgNormalizedPowerCSD(validNormCSD), alignedXCSD, 'k', 'LineWidth', 2); % Normalized power CSD
title('Normalized Aligned CSD Alpha Amp.');
xlabel('Normalized Average Power');
ylabel('Distance from inversion (mm)');
hold off;

% Plotting for LFP
subplot(1,2,2);
hold on;
fill([ciLowerNormLFP(validNormLFP); flipud(ciUpperNormLFP(validNormLFP))], ...
     [alignedXLFP, fliplr(alignedXLFP)], [0.9 0.9 0.9], 'EdgeColor', 'none');
plot(avgNormalizedPowerLFP(validNormLFP), alignedXLFP, 'k', 'LineWidth', 2); % Normalized power LFP
title('Normalized Aligned LFP Alpha Amp.');
xlabel('Normalized Average Power');
ylabel('Distance from inversion (mm)');
hold off;

%% Initialize array to store the channel indices with the highest power
highestChannelCSD = NaN(1, size(alignedPowerCSD, 2));
highestChannelLFP = NaN(1, size(alignedPowerLFP, 2));

% Find the maximum power channel index for each column
for col = 1:size(alignedPowerCSD, 2)
    [~, idxCSD] = max(alignedPowerCSD(:, col));
    [~, idxLFP] = max(alignedPowerLFP(:, col));
    highestChannelCSD(col) = idxCSD;
    highestChannelLFP(col) = idxLFP;
end

% Calculate distances from the reference channel
referenceIdx = totalRows / 2; % Assuming center is reference
channelDistances = ((1:totalRows) - referenceIdx) * 0.1; % Each channel spaced by 0.1 mm

% Determine non-zero data ranges for the y-axis
validRowsCSD = find(any(~isnan(alignedPowerCSD), 2));
validRowsLFP = find(any(~isnan(alignedPowerLFP), 2));

% Create histograms (channel counts) using actual channel numbers
channelCountCSD = histcounts(highestChannelCSD, 0.5:1:(totalRows+0.5));
channelCountLFP = histcounts(highestChannelLFP, 0.5:1:(totalRows+0.5));

%%
% Plotting histograms as horizontal bar graphs, with distances
figure;
subplot(1, 2, 1);
barh(channelCountCSD, 'FaceColor', 'b');
title('Max Power Channel (CSD)');
xlabel('Count');
ylabel('Distance from Inversion (mm)');
%set(gca, 'YDir', 'reverse', 'YTick', validRowsCSD, 'YTickLabel', flip(channelDistances(validRowsCSD)));  % Show distance from reference
%set(gca, 'YDir','YTick', validRowsCSD);
set(gca,'YTickLabel',(((20:5:45)-33)*0.1))
axis tight;
ylim([min(validRowsCSD)-1 max(validRowsCSD)+1]); % Adjust y-axis limits based on non-zero data

subplot(1, 2, 2);
barh(channelCountLFP, 'FaceColor', 'r');
title('Max Power Channel (LFP)');
xlabel('Count');
ylabel('Distance from Inversion (mm)');
%set(gca, 'YDir', 'reverse', 'YTick', validRowsLFP, 'YTickLabel', flip(channelDistances(validRowsLFP)));  % Show distance from reference
%set(gca, 'YDir','YTick', validRowsLFP);
set(gca,'YTickLabel',(((20:5:45)-33)*0.1))
axis tight;
ylim([min(validRowsLFP)-1 max(validRowsLFP)+1]); % Adjust y-axis limits based on non-zero data



%% Finalization
toc;