%%
clear data datavis plotcolor trig LFP CSD MUA maximum maximum_previous gimmeplots 
avgResp_aud_norm = avgResp_aud/max(avgResp_aud);
avgResp_vis_norm = AvgResps_vis/max(avgResp_aud);
%% Save the workspace
workspace_filename = strcat(paths{1,1}, '_workspace.mat');
save(workspace_filename,'-v7.3');

% Save the current figure
fig_filename = strcat(paths{1,1}, '_figure.fig');  % Modify as needed
saveas(gcf, fig_filename);