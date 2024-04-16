% SU analysis shell script 
close all;
clear;clc;
spikingDataFile = ('C:\Users\cmackey\Documents\AttnTuning\kk046\2-kk045046013@sortedv2.mat');
trigDataFile = ('C:\Users\cmackey\Documents\AttnTuning\kk046\2-kk045046013@trigs.mat');
channelIDs = [9,11,12];
unitIDs = [1:2];
preTime = .01; % second
postTime = 0.05; % seconds
sampleRate = 9000; % Hz
trigch = 1; % triggers 

% basic plots of SUA
plotSUData(spikingDataFile, trigDataFile, channelIDs, unitIDs, preTime, postTime, sampleRate,trigch)

% correlate SUA
[corrMatrix,fanoFactors,tickLabels] = calculateSpikeCountCorrelation(spikingDataFile, trigDataFile, channelIDs, unitIDs, preTime, postTime, sampleRate,trigch)

%% save stuff
% New section for adjusted file naming and saving
% Construct a string representation of the channel IDs for file naming
channelIDsStr = strjoin(arrayfun(@num2str, channelIDs, 'UniformOutput', false), '_');

% Extract base file name and folder path from spikingDataFile
[filepath, name, ~] = fileparts(spikingDataFile);
outputBaseFileName = sprintf('%s_Ch%s', name, channelIDsStr);

% Adjusting and saving figures with new file names
figureWidth = 1600; % Width in pixels
figureHeight = 900; % Height in pixels

% Save Raster Plot Figure
fig1 = figure(1);
set(fig1, 'Position', [100, 100, figureWidth, figureHeight]);
rasterFileName = fullfile(filepath, sprintf('%s_Rasters.fig', outputBaseFileName));
saveas(fig1, rasterFileName);

% Save Correlation Heatmap Figure
fig3 = figure(3);
correlationFileName = fullfile(filepath, sprintf('%s_Correlations.fig', outputBaseFileName));
saveas(fig3, correlationFileName);

% Save Correlation Info and Fano Factors to Text File

outputFileName = fullfile(filepath, sprintf('%s_correlationInfo.txt', outputBaseFileName));

% Open file for writing
fileID = fopen(outputFileName, 'w');

% Check if file was opened successfully
if fileID == -1
    error('Failed to open file for writing: %s', outputFileName);
end

% Write data to file
fprintf(fileID, 'Channel IDs: %s\n', mat2str(channelIDs));
fprintf(fileID, 'Unit IDs: %s\n', mat2str(unitIDs));
fprintf(fileID, 'Sample Rate: %d Hz\n', sampleRate);
fprintf(fileID, 'Trigger Channel: %d\n', trigch);
fprintf(fileID, 'Correlation Matrix:\n');

% Write the correlation matrix
for i = 1:size(corrMatrix, 1)
    fprintf(fileID, '%f\t', corrMatrix(i, :));
    fprintf(fileID, '\n');
end

% Add a separator or header to indicate the start of Fano Factors data
fprintf(fileID, '\nFano Factors for Each Channel-Unit Combination:\n');

% Write the Fano Factors
for i = 1:length(fanoFactors)
    % Assuming tickLabels already corresponds to the channel-unit combinations
    fprintf(fileID, '%s: %f\n', tickLabels{i}, fanoFactors(i));
end

% Close the file
fclose(fileID);
