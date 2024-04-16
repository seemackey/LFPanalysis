function troughTimes = plotSUData_findTroughs(spikeWaveforms)
    % This function finds the indices of the troughs in spike waveforms
    % Assuming the minimum voltage value indicates the trough
    [~, troughTimes] = min(spikeWaveforms, [], 2);
end