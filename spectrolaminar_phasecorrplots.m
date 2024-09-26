%% SPECTROLAMINAR PHASE CORRELATION PLOTS
% CHASE M 2024
%% Initialization
tic;
clear;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\VisCtxData\freev\cont\Alpha\'; % Specify directory
myfiles = dir(fullfile(mydir, '*AveragePowerPhaseData.mat'));  % Adjust the pattern to match your new data files

%% Process each file

channelsDesc = [];

for loopct = 1:length(myfiles)
    close all
    clear epochPhaseCSD epochPhaseLFP
    % Load the data file
    dataFilePath = fullfile(mydir, myfiles(loopct).name);
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    load(dataFilePath, 'epochPhaseCSD', 'epochPhaseLFP');  % Assuming phase data is stored as (channel, trial, time)

    % Extract the midpoint index for phase data
    midIndex = floor(size(epochPhaseCSD, 3) / 2) + 1;

    % Extract phase data at midpoint across all trials
    midPhaseCSD = [];
    midPhaseLFP = [];
    midPhaseCSD = cat(2, midPhaseCSD, squeeze(epochPhaseCSD(:, :, midIndex)));  % Concatenate horizontally across files
    midPhaseLFP = cat(2, midPhaseLFP, squeeze(epochPhaseLFP(:, :, midIndex)));


    %% Correlation calculations
    % Correlation within CSD channels
    corrMatrixCSD = corr(midPhaseCSD');

    % Correlation within LFP channels
    corrMatrixLFP = corr(midPhaseLFP');

    % Correlation between CSD and LFP channels
    corrMatrixCSDLFP = corr(midPhaseCSD', midPhaseLFP');

    
    
    %% Plotting correlation matrices
    
    % Define filenames for saving the figures and data
    % Create a new directory for figures if it doesn't exist
    figuresDir = fullfile(mydir, 'PhaseCorr');
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end

    % Update figure file paths to new directory
    figureFileNameCSDPhase = fullfile(figuresDir, [basefilename(1:end-4), '_PhaseCorr.fig']);


    dataFileName = fullfile(figuresDir, [basefilename(1:end-4), '_PhaseCorrData.mat']);
    
    % FIGGLES here
    figure;
    subplot(1, 3, 1);
    imagesc(corrMatrixCSD);
    title('Phase Corr. (CSD)');
    xlabel('Channels');
    ylabel('Channels');
    colorbar;
    axis square;
    caxis([-1 1]);
    

    subplot(1, 3, 2);
    imagesc(corrMatrixLFP);
    title('Phase Corr. (LFP)');
    xlabel('Channels');
    ylabel('Channels');
    colorbar;
    axis square;
    caxis([-1 1]);

    subplot(1, 3, 3);
    imagesc(corrMatrixCSDLFP);
    title('Phase Corr. (CSD to LFP)');
    xlabel('CSD Channels');
    ylabel('LFP Channels');
    colorbar;
    axis square;
    caxis([-1 1]);
    saveas(gcf, figureFileNameCSDPhase);

    %% Save the data
    save(dataFileName, 'corrMatrixCSD', 'corrMatrixLFP', 'corrMatrixCSDLFP','epochPhaseCSD', 'epochPhaseLFP');


end
%% Finalization
toc;
