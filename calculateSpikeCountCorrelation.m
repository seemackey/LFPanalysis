function [corrMatrix,fanoFactors,tickLabels] = calculateSpikeCountCorrelation(spikingDataFile, trigDataFile, channelIDs, unitIDs, preTime, postTime, sampleRate,trigch)
    % Load spiking data
    data = load(spikingDataFile);
    channelNames = fieldnames(data);
    channels = struct();
    for i = 1:length(channelNames)
        if startsWith(channelNames{i}, 'Channel')
            channels.(channelNames{i}) = data.(channelNames{i});
        end
    end

    % Load trigger data and filter by trigType == 1
    trigData = load(trigDataFile);
    trig = trigData.trig;
    [trigType, trigTimes] = EphysExtractTrigs(trig, trigch, sampleRate);
    filterIndex = trigType == 1;
    trigTimesFiltered = trigTimes(filterIndex);
    trigTimesSeconds = trigTimesFiltered / sampleRate;

    % Preallocate matrix to store spike counts
    numChannels = length(channelIDs);
    numUnits = length(unitIDs);
    numTriggers = length(trigTimesSeconds);
    spikeCounts = zeros(numChannels * numUnits, numTriggers);

    % Calculate spike counts for each channel and unit
    tickLabels = {}; % Initialize tickLabels to store channel-unit combination labels
    countIndex = 1;
    for c = channelIDs
        for u = unitIDs
            channelName = sprintf('Channel%02d', c);
            if isfield(channels, channelName)
                spikes = channels.(channelName);
                unitSpikes = spikes(spikes(:, 2) == u, :);
                spikeTimesSeconds = unitSpikes(:, 3); % Spike times
                
                for tIdx = 1:numTriggers
                    t = trigTimesSeconds(tIdx);
                    spikeCounts(countIndex, tIdx) = sum(spikeTimesSeconds >= t - preTime & spikeTimesSeconds <= t + postTime);
                end
                
                % Add label for this channel-unit combination
                tickLabels{end+1} = sprintf('C%02dU%d', c, u);
                countIndex = countIndex + 1;
            end
        end
    end

    % Spike Count Visualization
    figure;
    totalSpikeCounts = sum(spikeCounts, 2); % Sum across triggers for each channel-unit combination
    bar(totalSpikeCounts);
    title('Total Spike Counts for Each Channel-Unit Combination');
    xlabel('Channel-Unit Combination');
    ylabel('Total Spike Count');
    xticks(1:length(tickLabels));
    xticklabels(tickLabels);
    xtickangle(45);

    % Compute Pearson correlation matrix
    corrMatrix = corr(spikeCounts');

    % Plot heatmap of correlation matrix
    figure;
    subplot(1,2,1)
    imagesc(corrMatrix);
    title('Pearson Correlation of Spike Counts');
    xlabel('Channel-Unit Combination');
    ylabel('Channel-Unit Combination');
    colorbar;
    axis square;
    xticks(1:length(tickLabels));
    xticklabels(tickLabels);
    yticks(1:length(tickLabels));
    yticklabels(tickLabels);
    xtickangle(45);
    
    
    % Assuming spikeCounts is numUnits*numChannels x numTriggers
    fanoFactors = var(spikeCounts, 0, 2) ./ mean(spikeCounts, 2);
    
    % Assuming you already have a figure for the correlation matrix
    % Create a new subplot for Fano Factors
    subplot(1, 2, 2); % Adjust subplot grid as needed
    bar(fanoFactors);
    title('Fano Factor for Each Channel-Unit Combo');
    xlabel('Channel-Unit Combo');
    ylabel('Fano Factor');
    xticks(1:length(tickLabels));
    xticklabels(tickLabels);
    xtickangle(45);


end
