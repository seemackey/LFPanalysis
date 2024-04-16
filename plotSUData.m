function plotSUData(spikingDataFile, trigDataFile, channelIDs, unitIDs, preTime, postTime, sampleRate,trigch)
    % Load spiking data
    data = load(spikingDataFile);
    channelNames = fieldnames(data);
    channels = struct();
    for i = 1:length(channelNames)
        name = channelNames{i};
        if startsWith(name, 'Channel')
            channels.(name) = data.(name);
        end
    end
    
    % Load trigger data
    trigData = load(trigDataFile);
    trig = trigData.trig;
    [trigType, trigTimes] = EphysExtractTrigs(trig, trigch, sampleRate);
    % Filter for trigType == 1
    filterIndex = trigType == 1;
    trigTimesFiltered = trigTimes(filterIndex);
    trigTimesSeconds = trigTimesFiltered / sampleRate;
    
    
    % Figure setup
    numChannels = length(channelIDs);
    numUnits = length(unitIDs);
    plotsPerUnit = 3; % For Raster, PSTH, and Waveform
    figure;
    
    for cIdx = 1:numChannels
        c = channelIDs(cIdx);
        for uIdx = 1:numUnits
            u = unitIDs(uIdx);
            channelName = sprintf('Channel%02d', c);
            if isfield(channels, channelName)
                spikes = channels.(channelName);
                unitSpikes = spikes(spikes(:, 2) == u, :);
                spikeTimesSeconds = unitSpikes(:, 3); % Spike times
                
                % Raster Plot
                subplot(numChannels*numUnits, plotsPerUnit, (cIdx-1)*numUnits*plotsPerUnit + (uIdx-1)*plotsPerUnit + 1);
                hold on;
                for t = trigTimesSeconds
                    eventSpikes = spikeTimesSeconds(spikeTimesSeconds >= t - preTime & spikeTimesSeconds <= t + postTime) - t;
                    plot(eventSpikes, ones(size(eventSpikes)) * t, 'k.', 'MarkerSize', 5);
                end
                hold off;
                title(sprintf('Raster: Ch %02d, Unit %d', c, u));
                xlabel('Time (s)');
                ylabel('Trial');
                xlim([-preTime postTime]);
                
                % PSTH Plot
                subplot(numChannels*numUnits, plotsPerUnit, (cIdx-1)*numUnits*plotsPerUnit + (uIdx-1)*plotsPerUnit + 2);
                binSize = 0.002; % 
                edges = -preTime:binSize:postTime;
                psthData = zeros(length(edges)-1, 1);
                for t = trigTimesSeconds
                    eventSpikes = spikeTimesSeconds(spikeTimesSeconds >= t - preTime & spikeTimesSeconds <= t + postTime) - t;
                    counts = histcounts(eventSpikes, edges);
                    psthData = psthData + counts';
                end
                bar(edges(1:end-1), psthData / length(trigTimesSeconds) / binSize, 'FaceColor', 'k');
                title(sprintf('PSTH: Ch %02d, Unit %d', c, u));
                xlabel('Time (s)');
                ylabel('Rate (spikes/s)');
                xlim([-preTime postTime]);
                
                % Waveform Plot
                subplot(numChannels*numUnits, plotsPerUnit, (cIdx-1)*numUnits*plotsPerUnit + (uIdx-1)*plotsPerUnit + 3);
                waveforms = unitSpikes(:, 5:end);
                avgWaveform = mean(waveforms, 1);
                plot(avgWaveform, 'LineWidth', 2);
                title(sprintf('Waveform: Ch %02d, Unit %d', c, u));
                xlabel('Sample');
                ylabel('Amplitude (uV)');
                spikeWidth = mean(unitSpikes(:, 4)); % Spike width
                text(max(xlim)*0.8, max(ylim)*0.8, sprintf('Width: %.2f ms', spikeWidth), 'FontSize', 8);
            end
        end
    end
end