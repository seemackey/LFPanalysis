%% Loops through the ephys extract function and gets out LFP data

clear
close all

% make this 1 if you want to see acoustic responses 
gimmeplots=1;

% select trigger, 1 for aud, 3 for visual, 4 for sacc on, 5 for sacc off
trigch=1;

%% path info 

%  paths = {'D:\oldtono\contproc\1-gt030031021@os_eye06_30'
% 'D:\oldtono\contproc\1-gt030031021@os_eye06_30'};

%  paths = {
%      '/Volumes/14TB_USB_01/jitter/_missing/ma041/2-ma041042041@os.mat';
% };
 
  paths = {
'/Volumes/14TB_USB_01/jitter/ma025_01/1-ma025026043@os.mat';
'/Volumes/14TB_USB_01/jitter/ma025_01/1-ma025026044@os.mat';
'/Volumes/14TB_USB_01/jitter/ma025_01/1-ma025026045@os.mat';
'/Volumes/14TB_USB_01/jitter/ma025_01/1-ma025026046@os.mat';};
 
 %/Volumes/16TB_003/dyneyep/bbn/contproc/2-rb069070030@os_eye06_20.mat bbn

% mac path example
% paths = {'/Volumes/Samsung03/bbn/1-bu015016038@os.mat'
%     '/Volumes/Samsung03/bbn/2-bu015016038@os.mat'};


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\1-bu019020039@os_eye06_20.mat';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\2-bu019020039@os_eye06_20.mat'; };

%% loop through the function to extract data and find outliers
all=cell(length(paths));
tic
for loopct = 1:1:length(paths)
    
    [CSD, LFP, trig, trigtype, MUA,epoch_tframe] =  EphysExtractFxn(paths{loopct,1},trigch); % extract data and stim triggers

    %unclean data
    all{1,loopct} = num2cell(CSD); 
    all{2,loopct} = num2cell(LFP);
    all{3,loopct} = num2cell(MUA);
    all{4,loopct} = num2cell(trig); % trig times
    all{5,loopct} = num2cell(trigtype); % trigger types (e.g. tone freqs, LED types)
    
    [~, ~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig, MUA); % finding artifacts
    
    outliers{:,loopct} = outlieridxtmp; %at the end this has idx of outliers from both recordings

    
end
toc

%% remove outliers across two sites (only works for two right now)
outliers_mutual = [outliers{:,1};outliers{:,length(paths)}];
outliersunique = unique(outliers_mutual);

for outlierct = 1:length(paths)
  
    all{1,outlierct}(:,outliersunique,:) = []; 
    all{2,outlierct}(:,outliersunique,:) = []; 
    all{3,outlierct}(:,outliersunique,:) = []; 
    all{4,outlierct}(:,outliersunique,:) = []; 
    all{5,outlierct}(outliersunique,:) = [];
    
end
    
%% Plot CSD, MUA, tuning
chmua1=7;
chmua2=7;
% tpmua1=round(length(MUA(1,1,:))/2); %using halfway point because I've been epoching 250 ms before stim and 250 after
% tpmua2=length(MUA(1,1,:));
tpmua1=55; %using halfway point because I've been epoching 250 ms before stim and 250 after
tpmua2=90;


site1LFP = cell2mat(all{2,1}(:,:,:)); %using LFP for site 1 (mgb)
%site2CSD = cell2mat(all{1,length(paths)}(:,:,:)); % 2nd site, ctx
site2CSD = cell2mat(all{1,1}(:,:,:)); % 2nd site, ctx
MUA = cell2mat(all{3,1}(:,:,:)); % MUA from second site
trigtype = cell2mat(all{5,1}(:,:)); % 

%     for i=1:size(MUA,1)
%         for iii=1:size(MUA,2)
%             for times_smooth=1:2
%              MUA(i,iii,:)=smooth(squeeze(MUA(i,iii,:)),'moving',5);%smoothed data
%             end
%         end
%     end
%     
%     for i=1:size(site2CSD,1)
%         for iii=1:size(site2CSD,2)
%             for times_smooth=1:2
%              site2CSD(i,iii,:)=smooth(squeeze(site2CSD(i,iii,:)),'moving',5);%smoothed data
%             end
%         end
%     end
    
    
if gimmeplots == 1
    
    % check LFPs
%     figure
%     axpos=[0.1 0.1 0.8 0.8];
%     figureax1a=axes('Position',axpos);
%     [cax2] = csd_maker_no_subplot07(squeeze(mean(site1LFP(:,:,:),2)),(-30:2:250),1,[-20 250],[0 0],[],axpos,figureax1a);
%     colormap(figureax1a,flipud(jet))
%     figtitle = ['scale: ' num2str(LFPscale)];
%     title(figtitle)
    site2CSD = cell2mat(all{2,4}(:,:,:)); % 2nd site, ctx
    figure
    axpos=[0.1 0.1 0.8 0.8];
    figureax1a=axes('Position',axpos);
    [LFPscale] = csd_maker_no_subplot07(squeeze(mean(site2CSD(3:19,:,:),2)),epoch_tframe(1):1:epoch_tframe(2),1,[-50 300],[0 0],[-2.15 2.15],axpos,figureax1a);
    colormap(figureax1a,flipud(jet))
    figtitle = ['scale: ' num2str(LFPscale)];
    title(figtitle)
    % COMMENT IN FOR DEPTH Y AXIS
    
    tickdepth = 0:0.1:2.1;
%     labels = {'0','','0.2','','0.4','','0.6','','0.8','',...
%        '1.0','','1.2','','1.4','','1.6','','1.8','','2.0','',};
       labels = [];
    yticklabels(labels)
    %ylabel('Depth (mm)')
    % COMMENT IN FOR DEPTH Y AXIS
    
    xlabel('Time (ms)')
    
% MUA MUA MUA MUA MUA MUA MUA
%     figure
%     axpos=[0.1 0.1 0.8 0.8];
%     figureax1a=axes('Position',axpos);
%     [MUAscale] = csd_maker_no_subplot07(squeeze(mean(MUA(:,:,:),2)),epoch_tframe(1):1:epoch_tframe(2),0,[-50 300],[0 0],[],axpos,figureax1a);
%     colormap(figureax1a,hot)
%     fig2title = ['scale: ' num2str(MUAscale)];
%     title(fig2title)
    
    % COMMENT IN FOR DEPTH Y AXIS
    tickdepth = 0:0.1:2.1;
%     labels = {'0','','0.2','','0.4','','0.6','','0.8','',...
%        '1.0','','1.2','','1.4','','1.6','','1.8','','2.0','',};
       labels = [];
    yticklabels(labels)
    %ylabel('Depth (mm)')
    % COMMENT IN FOR DEPTH Y AXIS
    
    xlabel('Time (ms)')
    
    
%%     % now check tuning
%     MUAavg = mean(mean(MUA(chmua1:chmua2,:,tpmua1:tpmua2),1),3)'; % average MUA across chs and whole epoch
%     uniquetrigs = unique(trigtype);
%     
%     
%     % sort by trigger type (e.g. tone frequency or led type
%     for trig_ct = 1:1:length(uniquetrigs)
%         
%         currtrig = trigtype(:,1)==uniquetrigs(trig_ct,1);
%         MUAavgsort(trig_ct,1) = mean(MUAavg(currtrig==1,1));
%         MUAstdsort(trig_ct,1) = std(MUAavg(currtrig==1,1));
%         
%     end
%     freqxaxis = [0.35,0.5,0.7,1,1.4,2,2.8,4,5.6,8,11.3,16,22.6,32];
%     figure
%     subplot(1,2,1)
%     errorbar(uniquetrigs,MUAavgsort(:,1),MUAstdsort)
%     xlabel('Trigger Type')
%     title('Tuning Curve')
%     %ylim([-1,3])
%     if length(uniquetrigs)>2
%         subplot(1,2,2)
%         semilogx(freqxaxis,(MUAavgsort((1:end-1),1)/max(MUAavgsort((1:end-1),1))))
%         xlabel('Tone Freq (kHz)')
%         title('Normalized tuning curve')
%     else
%         subplot(1,2,2)
%         plot(uniquetrigs,(MUAavgsort(:,1)/max(MUAavgsort(:,1))))
%         xlabel('Tone Freq (kHz)')
%         title('Normalized tuning curve')
%     end
    
end

%% calculate confidence intervals of MUA and significance of response
% numchans = length(all{1,1}(:,1,1));
% numtrs = length(all{1,1}(1,:,1));
% numtps = length(squeeze(all{1,1}(1,1,:)));
% ci_stp = 1; % step size, multiply by 2 to get milliseconds
% numboots = 500; % number of resamples
% tpmua2=800;
% % loop through channels
% for ci_ch_ct = 1:1:numchans
%     for ci_tp = 1:ci_stp:1001 % loop through time points
%         MUAtmp = squeeze(MUA(ci_ch_ct,:,ci_tp)); % get the data from a time point
%         ci_tmp = bootci(numboots,@mean,MUAtmp'); % get bootstrap derived CIs
%         cimean(ci_ch_ct,ci_tp) = mean(MUAtmp); % record mean
%         cilo(ci_ch_ct,ci_tp) = ci_tmp(1,1); % lower CI
%         cihi(ci_ch_ct,ci_tp) = ci_tmp(2,1); % upper CI
%         
%         if cilo(ci_ch_ct,ci_tp)>0
%             sigtimes(ci_ch_ct,ci_tp) = 1;
%         elseif cihi(ci_ch_ct,ci_tp)<0
%             sigtimes(ci_ch_ct,ci_tp) = -1;
%         else
%             sigtimes(ci_ch_ct,ci_tp) = 0;
%         end
%     end
% end
% 
% figure
% plotstp = ci_stp*1;
% tpmua2plot=1001;
% for ciplotct = 1:1:numchans
%     subplot(numchans,1,ciplotct)
%     plot(-200:plotstp:tpmua2,cihi(ciplotct,1:ci_stp:tpmua2plot))
%     hold on
%     plot(-200:plotstp:tpmua2,cilo(ciplotct,1:ci_stp:tpmua2plot))
%     hold on
%     plot(-200:plotstp:tpmua2,cimean(ciplotct,1:ci_stp:tpmua2plot))
%     ylim([-0.75,0.75])
% end
% 
% figure
% plotstp = ci_stp*1;
% increases=[0:0.9:numchans*0.9];
% for ciplotct = 1:1:numchans
%     
%     plot(-200:plotstp:tpmua2,cihi(ciplotct,1:ci_stp:tpmua2plot)+increases(ciplotct),'r')
%     hold on
%     plot(-200:plotstp:tpmua2,cilo(ciplotct,1:ci_stp:tpmua2plot)+increases(ciplotct),'r')
%     hold on
%     plot(-200:plotstp:tpmua2,cimean(ciplotct,1:ci_stp:tpmua2plot)+increases(ciplotct),'b')
%     hold on
%     plot(-200:plotstp:tpmua2,(zeros(length(-200:plotstp:tpmua2))+increases(ciplotct)),'k')
%     
% end
% ylim([-0.5,20])

%% %% plot CSD across trials
% numchans = length(all{1,1}(:,1,1));
% numtrs = length(all{1,1}(1,:,1));
% numtps = length(squeeze(all{1,1}(1,1,:)));
% ci_stp = 1; % step size, multiply by 2 to get milliseconds
% numboots = 500; % number of resamples
% chan = 10;
% 
% % loop through channels
% for ci_tr_ct = 1:1:numtrs
%     
%         CSDforplot(ci_tr_ct,:) = squeeze(site2CSD(chan,ci_tr_ct,:)); % get the data from a time point
% 
% end
% 
% figure
% plotstp = ci_stp*2;
% for ciplotct = 10:1:20
%     subplot(11,1,ciplotct-9)
% 
%     hold on
%     plot(-250:2:251,CSDforplot(ciplotct,:))
%      ylim([-20,20])
% end
% 
% 
