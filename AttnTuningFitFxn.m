function [maximum,gaussianFit,freqs] = AttnTuningFitFxn(all,selectedChannel,selectedTimeWindow,maximum_previous,subplotArgs,plotcolor)
%% inputs %%
% all, which contains all data from the extract csd loop %
% for now, "all" needs to be organized in a certain order specified in the
% extract CSD script %
%% outputs %% are just figures with all important data in them %
% chase mackey 2024 %

    %% Define channel, time window, and trigger type of interest

    triggerTypeOfInterest = 1;  % Change this to the desired trigger type

    % Initialize variables for storing average and CI
    averageResponses = [];
    ciResponses = [];

    % Loop through each file in "all"
    for fileIdx = 1:size(all, 2)
        % Extract MUA for the selected channel, time window, and trigger type
        selectedTrials = find(cell2mat(all{5, fileIdx}) == triggerTypeOfInterest);
        selectedMUA = cell2mat(all{3, fileIdx}(selectedChannel, selectedTrials, selectedTimeWindow));

        % Reshape the data for bootstrapping
        reshapedMUA = reshape(selectedMUA, size(selectedMUA, 1), []);

        % Perform bootstrapping and calculate confidence intervals
        bootstrappedMeans = bootstrp(1000, @mean, reshapedMUA');
        ci = prctile(bootstrappedMeans, [10, 90]);

        % Calculate average across time and channels
        avg = mean(selectedMUA(:));  % Average across both time and channels

        % Store results
        averageResponses = [averageResponses, avg];
        ciResponses = [ciResponses, ci];
    end

    % Plot average responses and error bars for each file
    %colors = parula(size(all, 2));
    hold on
    subplot(subplotArgs(1), subplotArgs(2), subplotArgs(3));
    for fileIdx = 1:size(all, 2)
        % Extract the corresponding frequency for each file
        frequency = cell2mat(all{6, fileIdx});
        freqs(fileIdx) = cell2mat(all{6, fileIdx});

        % Plot the average response and error bars
        hold on
        errorbar(frequency, averageResponses(fileIdx), ciResponses(1, fileIdx) - averageResponses(fileIdx), ciResponses(2, fileIdx) - averageResponses(fileIdx), 'o','Color', plotcolor);
    end
    xlabel('Frequency');
    ylabel('Average MUA');
    title('Average MUA Responses');
    set(gca, 'XScale', 'log');  % Set x-axis to log scale

    % Check if there are at least 4 different frequencies
    uniqueFrequencies = unique(cell2mat(cellfun(@cell2mat, all(6, :), 'UniformOutput', false)));

    if numel(uniqueFrequencies) >= 4
        if maximum_previous == 0
            % Normalize responses to the highest frequency
            normalizedResponses = averageResponses / max(averageResponses);
            maximum = max(averageResponses); % saving this as an ouput of the fxn
        else
            normalizedResponses = averageResponses / maximum_previous; % norm to previous data's max
            maximum = 0; % we didn't use maximum, we used max prev.
        end
        
        % Fit a Gaussian function
        gaussianFit = fit(freqs', normalizedResponses', 'gauss1');

        % Display fit parameters
        fprintf('Amplitude: %.3f\n', gaussianFit.a1);
        fprintf('Mean: %.3f\n', gaussianFit.b1);
        fprintf('Width: %.3f\n', gaussianFit.c1 / sqrt(2));
        fprintf('Asymptote: %.3f\n', gaussianFit.a1 + gaussianFit.c1 / sqrt(2));

      
        % Plot the normalized responses and Gaussian fit
        hold on
        % Use subplotArgs to specify the subplot configuration
        subplot(subplotArgs(1), subplotArgs(2), subplotArgs(3)+1);
%         plot(freqs, normalizedResponses, 'o', 'Color', plotcolor);
%         hold on;
        plot(gaussianFit, freqs, normalizedResponses);
        h = get(gca, 'Children');
        set(h(1:2), 'Color', plotcolor);
        set(gca, 'XScale', 'log');  % Set x-axis to log scale
        xlabel('Frequency');
        ylabel('Normalized Average MUA');
        title('Normalized Average MUA Responses with Gaussian Fit');
        legend('Location', 'Best');
        % Annotate fit parameters on the figure
        %annotation('textbox', 'Position', [0.1, 0.8, 0.1, 0.1], 'String', sprintf('Amplitude: %.3f\nMean: %.3f\nWidth: %.3f\nAsymptote: %.3f', gaussianFit.a1, gaussianFit.b1, gaussianFit.c1 / sqrt(2), gaussianFit.a1 + gaussianFit.c1 / sqrt(2)), 'FitBoxToText', 'on', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    else
        disp('< 4 freqs, not enough to curve fit')
    end

end