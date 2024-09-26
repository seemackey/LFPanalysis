%% spectrolaminar, looking at crossing point of two power trends

%% Initialization
clc;
clear;
close all;

% Directory settings
datadir1 = 'E:\spectrolaminar\AttnData\core\cont\Alpha\';
datadir2 = 'E:\spectrolaminar\AttnData\core\cont\FiguresGamma\';

% Get list of data files from both directories
dataFiles1 = dir(fullfile(datadir1, '*AveragePowerPhaseData.mat'));
dataFiles2 = dir(fullfile(datadir2, '*AveragePowerPhaseData.mat'));

%% Prepare to process files
% Convert directory 2 filenames to a list for easier comparison
dataFiles2Names = {dataFiles2.name};

%% Processing files and extracting data
for k = 1:length(dataFiles1)
    fileName1 = dataFiles1(k).name;

    % Find matching file in the second directory
    idx = find(strcmp(dataFiles2Names, fileName1));
    if isempty(idx)
        continue; % Skip if no matching file is found
    end

    % Load data from both directories for the matching file
    data1 = load(fullfile(datadir1, fileName1), 'avgPowerLFPDesc');
    data2 = load(fullfile(datadir2, dataFiles2Names{idx}), 'avgPowerLFPDesc');

    % Normalizing LFP power data to the maximum channel in each file
    normPowerLFP1 = data1.avgPowerLFPDesc / max(data1.avgPowerLFPDesc);
    normPowerLFP2 = data2.avgPowerLFPDesc / max(data2.avgPowerLFPDesc);


    
        %% Channels for plotting (assuming all data files have the same number of channels)
    channels = 1:length(normPowerLFP1);
    both = [normPowerLFP1,normPowerLFP2];
    
    pow_min = min(min(both));
    % Plotting the normalized LFP power trends for matching files
    figure; % Create a new figure for each matching pair
    hold on;
    plot(flip(normPowerLFP1), channels, 'b', 'LineWidth', 2); % Normalized power from the first directory
    plot(flip(normPowerLFP2), channels, 'r', 'LineWidth', 2); % Normalized power from the second directory
    legend('A-B', 'Gamma');
    xlabel('Normalized Power');
    ylabel('Channel');
    set(gca,'Xlim',[pow_min 1]);
    title(['Alpha-Gamma Cross: ', fileName1(1:13)]);
    set(gca, 'YDir', 'reverse'); % Reverse y-axis to have channel 1 at the top
    grid on;
    hold off;


    
end
