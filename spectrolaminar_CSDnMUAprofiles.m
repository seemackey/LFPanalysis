% make spatiotemporal profiles of a dir of files
% chase m 2024

%% Initialization
tic;
clear;
clearvars;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\VisCtxData\freev\cont\'; % Specify directory
figuresDir = fullfile(mydir, 'CSDnMUAprofiles'); % dir where the figs go
allFlipResults = struct();
pickchannels = 0; % do you want to manually select channels based on the CSD? 
loadchannels = 0; % do you want to load previously selected channels?
peterdataflag = 0; % peter and noelle's A1 and belt Data
trigch = 1; % "1" = trigger code of aud stim in peter data 

if loadchannels == 1
    % dir with flip results and channel selections
    SelChansDir = fullfile(mydir, 'Crosses_taper_FLIP'); 
end

% handle file extension differences
if peterdataflag == 1
    myfiles = dir(fullfile(mydir, '*@os*.mat')); % Get all aud files in struct
else
    myfiles = dir(fullfile(mydir, '*@oe*.mat')); % annie's vis data
end

if isempty(myfiles)
    disp('No files found!');
end

% Load AllFlipResults if it exists because we want to use our saved
% "Selected Channels variable" for the visual cortex data only!
if loadchannels == 1
    % resultsFile = fullfile(mydir, 'allFlipResults.mat');
    % if exist(resultsFile, 'file')
    %     load(resultsFile, 'allFlipResults');
    % else
    %     allFlipResults = struct();
    % end
    load('E:\spectrolaminar\VisCtxData\cont\FPnMUA\allFlipResults.mat');
    allFlipResults = aNewSelectedChannels;
end

for loopct = 1:length(myfiles)
    close all
    clear selectedChannels
    %% Data3D with dimensions channels, trials, epoched time
    % Load the LFP file
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    
    % Determine MUA file name
    muaFilename = strrep(basefilename, '@oe', '@om');  % Replace '@oe' with '@om' to find the MUA file
    muaFullFilename = fullfile(mydir, muaFilename);

    % Olive recordings had 200 micron spacing
    if peterdataflag == 0
        if basefilename(3) == 'O'
            interelectrodespacing = 0.2;
        else
            interelectrodespacing = 0.1;
        end
    else
        interelectrodespacing = 0.1;
    end
    
    %% Load the LFP data
    if peterdataflag == 1
        [~, ~, ~, ~, ~, ~, ~, lfp] = EphysExtractFxn(fullfilename, trigch);
    else
        load(fullfilename);
        newadrate = 1000;
        epoch_tframe = [-300 300]; 
        trig0 = anatrig{1};
        trig01 = [];
        for trigredxct = 1:length(trig0)
            trig01(trigredxct) = round(trig0(trigredxct) / (adrate / newadrate));
        end
        [lfp] = EphysEpochFxnSimple(cnt, trig01, epoch_tframe, newadrate, adrate);
    end
    
    %% Check if selected channels exist in allFlipResults
    if loadchannels == 1
        if isfield(allFlipResults, 'filename') && any(strcmp({allFlipResults.filename}, basefilename))
            currentFileIndex = find(strcmp({allFlipResults.filename}, basefilename));
            if isfield(allFlipResults(currentFileIndex).results, 'selectedChannels')
                selectedChannels = allFlipResults(currentFileIndex).results.selectedChannels;
            else
                selectedChannels = 1:1:size(lfp,1);
            end
        else
            selectedChannels = 1:1:size(lfp,1);  % Default to all channels if none selected
        end
    end

    %% Apply Selected Channels
    if ~exist('selectedChannels','var')
        selectedChannels = 1:size(lfp,1);
        lfp = lfp(selectedChannels, :, :);
    else
        lfp = lfp(selectedChannels,:,:);
    end
    %% Load MUA Data if it exists
    muaExists = isfile(muaFullFilename);
    if muaExists
        if peterdataflag == 1
            [~, ~, ~, ~, ~, ~, mua] = EphysExtractFxn(muaFullFilename, trigch);
        else
            load(muaFullFilename);
            [mua] = EphysEpochFxnSimple(cnt, trig01, epoch_tframe, newadrate, adrate);
        end
        % Apply the same selected channels to MUA
        mua = mua(selectedChannels, :, :);
    else
        mua = [];  % No MUA data available
    end

    %% CSD Calculation
    clear CSD
    CSD(:, :, :) = -diff(lfp, 2, 1); % CSD
    lfp = lfp(2:end-1, :, :);  % Adjust LFP to match the size after CSD calculation

    % If MUA exists, trim it to match the size after CSD calculation
    if ~isempty(mua)
        mua = mua(2:end-1, :, :);
    end

    % artifact reject
    [mua, outliermua] = MTF_rejectartifacts(mua,'median',2);
    [CSD, ~] = MTF_rejectartifacts(CSD,'median',2);
    [lfp, ~] = MTF_rejectartifacts(lfp,'median',2);


    baselineidx = 250:300;

% Initialize baseline-corrected arrays
lfp_bsl = [];
CSD_bsl = [];
mua_bsl = [];

% Loop through trials of LFP
for trct = 1:size(lfp, 2)
    for chct = 1:size(lfp, 1)
        lfp_bsl(chct, trct, :) = squeeze(lfp(chct, trct, :)) - squeeze(mean(lfp(chct, trct, baselineidx)));
    end
end

% Loop through trials of CSD
for trct = 1:size(CSD, 2)
    for chct = 1:size(CSD, 1)
        CSD_bsl(chct, trct, :) = squeeze(CSD(chct, trct, :)) - squeeze(mean(CSD(chct, trct, baselineidx)));
    end
end

% Loop through trials of MUA if it exists
if ~isempty(mua)
    for trct = 1:size(mua, 2)
        for chct = 1:size(mua, 1)
            mua_bsl(chct, trct, :) = squeeze(mua(chct, trct, :)) - squeeze(mean(mua(chct, trct, baselineidx)));
        end
    end
end


    %% plots
    timeVector = 290:500;
    numChannels = 1:size(lfp_bsl, 1);
    numChannelscsd = 1:size(CSD_bsl, 1);

    % Determine Color Axis for CSD
    csd_min = min(min(squeeze(mean(CSD_bsl(:, :, timeVector), 2))));
    csd_max = max(max(squeeze(mean(CSD_bsl(:, :, timeVector), 2))));
    csd_caxis = [-max(abs([csd_min, csd_max])) * 0.6, max(abs([csd_min, csd_max])) * 0.6];

    % Determine Color Axis for MUA (if available)
    if ~isempty(mua_bsl)
        mua_min = min(min(squeeze(mean(mua_bsl(:, :, timeVector), 2))));
        mua_max = max(max(squeeze(mean(mua_bsl(:, :, timeVector), 2))));
        mua_caxis = [-max(abs([mua_min, mua_max])) * 0.6, max(abs([mua_min, mua_max])) * 0.6];
    end

    figure('Position', [200, 200, 1200, 700]);

    % Plot CSD
    subplot(1, 3, 1);
    imagesc(timeVector, numChannelscsd, squeeze(mean(CSD_bsl(:, :, timeVector), 2)));
    title('Trial Avg. CSD');
    xlabel('Time (s)');
    ylabel('Channel');
    colormap(flipud(jet)); 
    colorbar;
    caxis(csd_caxis);

    % Plot MUA if available
    if ~isempty(mua_bsl)
        subplot(1, 3, 3);
        imagesc(timeVector, numChannelscsd, squeeze(mean(mua_bsl(:, :, timeVector), 2)));
        title('Trial Avg. MUA');
        xlabel('Time (s)');
        ylabel('Channel');
        colorbar;
        caxis(mua_caxis);
    end

    % Call the FindSigChannel function
    
    [sigchans,sigtimes] = FindSigChannel(CSD_bsl(:,:,timeVector));

    % Plot significant channels
    

    hold on;
    subplot(1,3,2)
    imagesc(timeVector, 1:size(CSD_bsl, 1), squeeze(mean(CSD_bsl(:,:,timeVector), 2)));
    title('Significant CSD Channels');
    xlabel('Time (ms)');
    ylabel('Channel');
    set(gca, 'YDir', 'reverse');
    colorbar;
    caxis(csd_caxis); % Adjust color axis for better visualization
    % Overlay significant activation times with circles
    hold on;
    [sig_ch, sig_time] = find(sigtimes < 0); % Find significant time points
    plot(timeVector(sig_time), sig_ch, 'o', 'MarkerSize', 4, 'Color', [1 0 0]); % Red circles for significant points
    hold off;

    %% Save the figure
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end

    figureFileName = fullfile(figuresDir, [basefilename(1:end-4), '_Profiles.fig']);
    saveas(gcf, figureFileName);
    saveas(gcf, [figureFileName(1:end-4), '.jpg']);
    close(gcf);

    disp(['Processed and saved: ', basefilename]);
end

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

toc;

function [cleaned_data, outlier_idx_unique] = MTF_rejectartifacts(data, method, threshold)
    % Artifact rejection function
    % data: Input data (channels x trials x time)
    % method: Method to identify outliers ('median' or 'mean')
    % threshold: Threshold for identifying outliers

    if nargin < 2
        method = 'median';  % Default method
    end
    
    if nargin < 3
        threshold = 3;  % Default threshold (z-score)
    end

    num_trials = size(data, 2);

    max_vals = zeros(1, num_trials);

    % Calculate the maximum mean absolute value for each trial
    for trial = 1:num_trials
        max_vals(trial) = max(mean(abs(squeeze(data(:, trial, :))), 1));
    end

    % Identify outliers based on the chosen method
    switch method
        case 'median'
            outlier_idx = abs(max_vals - median(max_vals)) > threshold * std(max_vals);
        case 'mean'
            outlier_idx = abs(max_vals - mean(max_vals)) > threshold * std(max_vals);
        otherwise
            error('Unsupported method. Use "median" or "mean".');
    end

    outlier_idx_unique = find(outlier_idx);

    % Remove outliers
    cleaned_data = data;
    cleaned_data(:, outlier_idx_unique, :) = [];

end

function [sigchans, sigtimes] = FindSigChannel(CSD)

    % Smooth the CSD data with a 5-sample moving average window
    smoothingWindow = 5;
    for ch = 1:size(CSD, 1)
        for tr = 1:size(CSD, 2)
            CSD(ch, tr, :) = smooth(squeeze(CSD(ch, tr, :)), smoothingWindow);
        end
    end

    % Calculates significance of response in CSD data
    % CSD is the CSD data (channels x trials x time points)
    
    numchans = size(CSD, 1);
    numtrs = size(CSD, 2);
    numtps = size(CSD, 3);
    ci_stp = 1; % Step size
    numboots = 500; % Number of resamples
    baseline_tp = 1:20; % Baseline time points (first 20)
    
    sigchannels = {}; % Initialize sigchannels
    sigtimes = zeros(numchans, numtps); % Initialize sigtimes

    % Calculate baseline CIs
    baseline_cilo = zeros(numchans, 1);
    baseline_cihi = zeros(numchans, 1);

    for ci_ch_ct = 1:numchans
        baseline_data = squeeze(CSD(ci_ch_ct, :, baseline_tp)); % Get baseline data
        baseline_ci_tmp = bootci(numboots, @mean, baseline_data(:)); % Get baseline CIs
        baseline_cilo(ci_ch_ct) = baseline_ci_tmp(1); % Baseline lower CI
        baseline_cihi(ci_ch_ct) = baseline_ci_tmp(2); % Baseline upper CI
    end

    % Loop through channels and time points
    for ci_ch_ct = 1:numchans
        new = 1;
        for ci_tp = 1:ci_stp:numtps % Loop through time points
            MUAtmp = squeeze(CSD(ci_ch_ct, :, ci_tp)); % Get the data from a time point
            ci_tmp = bootci(numboots, @mean, MUAtmp'); % Get bootstrap derived CIs
            cimean(ci_ch_ct, ci_tp) = mean(MUAtmp); % Record mean
            cilo(ci_ch_ct, ci_tp) = ci_tmp(1); % Lower CI
            cihi(ci_ch_ct, ci_tp) = ci_tmp(2); % Upper CI

            % Determine if the CSD signal is significant based on stricter criteria
            if cilo(ci_ch_ct, ci_tp) > baseline_cihi(ci_ch_ct) % Signal is above baseline upper CI
                sigtimes(ci_ch_ct, ci_tp) = 1;
            elseif cihi(ci_ch_ct, ci_tp) < baseline_cilo(ci_ch_ct) % Signal is below baseline lower CI
                sigtimes(ci_ch_ct, ci_tp) = -1;
            else
                sigtimes(ci_ch_ct, ci_tp) = 0;
            end
            
            % Identify significant channels
            if ci_tp > numtps / 2
                if all(sigtimes(ci_ch_ct, ci_tp - 9:ci_tp) < 0)
                    sigchannels{ci_ch_ct, 1} = ci_ch_ct;
                    new = 0;
                end
            end


        end
    end
    
    % Pull out the channels with long enough activation times
    sigchans = unique(cell2mat(sigchannels));
end

