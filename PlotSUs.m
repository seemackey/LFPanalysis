%% import plexon sorted unit data and plot psth

% Sample data: spike times and event trigger times
spikeTimes = Channel11(:,3);
eventTimes = trig_downsamp';

% Set parameters
binWidth = 1;  % Width of each time bin (in seconds)
timeWindow = 5;  % Total time window to consider (in seconds), before and after the event

% Initialize arrays to store spike counts for each time bin
edges = -timeWindow:binWidth:timeWindow;
spikeCountsBefore = zeros(1, length(edges) - 1);
spikeCountsAfter = zeros(1, length(edges) - 1);

% Loop through each event
for i = 1:length(eventTimes)
    % Find spikes within the time window before the event
    spikesBefore = spikeTimes(spikeTimes >= (eventTimes(i) - timeWindow) & spikeTimes < eventTimes(i));
    
    % Find spikes within the time window after the event
    spikesAfter = spikeTimes(spikeTimes >= eventTimes(i) & spikeTimes < (eventTimes(i) + timeWindow));
    
    % Compute spike counts in each time bin using a loop
    for j = 1:length(edges) - 1
        spikeCountsBefore(j) = spikeCountsBefore(j) + sum(spikesBefore >= edges(j) & spikesBefore < edges(j + 1));
        spikeCountsAfter(j) = spikeCountsAfter(j) + sum(spikesAfter >= edges(j) & spikesAfter < edges(j + 1));
    end
end

% Plot the peri-event time histogram
timeBins = linspace(-timeWindow, timeWindow, length(spikeCountsBefore) + 1);
figure;
subplot(2,1,1);
bar(timeBins(1:end-1), spikeCountsBefore, 'b');
title('Spike Rate Before Event');
xlabel('Time (s)');
ylabel('Spike Count');

subplot(2,1,2);
bar(timeBins(1:end-1), spikeCountsAfter, 'r');
title('Spike Rate After Event');
xlabel('Time (s)');
ylabel('Spike Count');


