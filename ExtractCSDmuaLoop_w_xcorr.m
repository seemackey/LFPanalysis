%% Loops through the ephys extract function and gets out LFP data

clear
close all

% make this 1 if you want to see acoustic responses 
gimmeplots=1;

% select trigger, 1 for aud, 3 for visual, 4 for sacc on, 5 for sacc off
trigch=3;

%% path info 

%  paths = {'D:\oldtono\contproc\1-gt030031021@os_eye06_30'
% 'D:\oldtono\contproc\1-gt030031021@os_eye06_30'};

 paths = {'D:\led\contproc\New folder\1-bu009010028@os_eye06_20.mat'
'D:\led\contproc\New Folder\2-bu009010028@os_eye06_20'};


% mac path example
% paths = {'/Volumes/Samsung03/bbn/1-bu015016038@os.mat'
%     '/Volumes/Samsung03/bbn/2-bu015016038@os.mat'};

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\1-bu019020037@os_eye06_20';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\2-bu019020037@os_eye06_20'};


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu035036\1-bu035036022@os_eye06_20';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu035036\2-bu035036022@os_eye06_20'};

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu001002\1-bu001002020@os.mat';
%    };

%tono
% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\oldtono\1-bu035036019@os';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\oldtono\2-bu035036019@os'};


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu035036\1-bu035036021@os_eye06_20';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu035036\2-bu035036021@os_eye06_20'};

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\oldtono\1-bu019020027@os';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\oldtono\2-bu019020027@os'};


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu009010\1-bu009010034@os';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu009010\2-bu009010034@os'};

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\chase\eye\1-bu00901002_eye06_20'
%     '\\NKI-LAKATOSLAB\lakatoslab_alpha\chase\eye\2-bu00901002_eye06_20'};

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\chase\eye\1-bu015016034@os_eye06_20';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\chase\eye\2-bu015016034@os_eye06_20'};


% paths = {'C:\Users\cmackey\Documents\MATLAB\1-bu009010034@os';
%     'C:\Users\cmackey\Documents\MATLAB\2-bu009010034@os'
%     };

% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\1-bu019020039@os_eye06_20.mat';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu019020\2-bu019020039@os_eye06_20.mat'; };

%% loop through the function to extract data and find outliers
all=cell(length(paths));
tic
for loopct = 1:1:length(paths)
    
    [CSD, LFP, trig, trigtype, MUA] =  EphysExtractFxn(paths{loopct,1},trigch); % extract data and stim triggers

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
chmua1=9;
chmua2=9;
tpmua1=round(length(MUA(1,1,:))/2); %using halfway point because I've been epoching 250 ms before stim and 250 after
tpmua2=length(MUA(1,1,:));


site1LFP = cell2mat(all{2,1}(:,:,:)); %using LFP for site 1 (mgb)
site2CSD = cell2mat(all{1,length(paths)}(:,:,:)); % 2nd site, ctx
MUA = cell2mat(all{3,1}(:,:,:)); % MUA from second site
trigtype = cell2mat(all{5,1}(:,:)); % 

if gimmeplots == 1
    
    % check LFPs
    figure
    axpos=[0.1 0.1 0.8 0.8];
    figureax1a=axes('Position',axpos);
    [cax2] = csd_maker_no_subplot07(squeeze(mean(site1LFP(:,:,:),2)),(-250:2:250),1,[-50 250],[0 0],[],axpos,figureax1a);
    colormap(figureax1a,flipud(jet))
    
    if length(paths)>1
        figure
        axpos=[0.1 0.1 0.8 0.8];
        figureax1a=axes('Position',axpos);
        [cax2] = csd_maker_no_subplot07(squeeze(mean(site2CSD(:,:,:),2)),(-250:2:250),1,[-50 250],[0 0],[],axpos,figureax1a);
        colormap(figureax1a,flipud(jet))
    end
    
figure
axpos=[0.1 0.1 0.8 0.8];
figureax1a=axes('Position',axpos);
[cax2] = csd_maker_no_subplot07(squeeze(mean(MUA(:,:,:),2)),(-250:2:250),0,[-50 250],[0 0],[],axpos,figureax1a);
colormap(figureax1a,hot)
    
    % now check tuning
    MUAavg = mean(mean(MUA(chmua1:chmua2,:,tpmua1:tpmua2),1),3)'; % average MUA across chs and whole epoch
    uniquetrigs = unique(trigtype);
    
    
    % sort by trigger type (e.g. tone frequency or led type
    for trig_ct = 1:1:length(uniquetrigs)
        
        currtrig = trigtype(:,1)==uniquetrigs(trig_ct,1);
        MUAavgsort(trig_ct,1) = mean(MUAavg(currtrig==1,1));
        MUAstdsort(trig_ct,1) = std(MUAavg(currtrig==1,1));
        
    end
    freqxaxis = [0.35,0.5,0.7,1,1.4,2,2.8,4,5.6,8,11.3,16,22.6,32];
    figure
    subplot(1,2,1)
    errorbar(uniquetrigs,MUAavgsort(:,1),MUAstdsort)
    xlabel('Trigger Type')
    title('Tuning Curve')
    if length(uniquetrigs)>2
        subplot(1,2,2)
        semilogx(freqxaxis,(MUAavgsort((1:end-1),1)/max(MUAavgsort((1:end-1),1))))
        xlabel('Tone Freq (kHz)')
        title('Normalized tuning curve')
    else
        subplot(1,2,2)
        plot(uniquetrigs,(MUAavgsort(:,1)/max(MUAavgsort(:,1))))
        xlabel('Tone Freq (kHz)')
        title('Normalized tuning curve')
    end
    
end

%% calculate confidence intervals of MUA and significance of response
numchans = length(all{1,1}(:,1,1));
numtrs = length(all{1,1}(1,:,1));
numtps = length(squeeze(all{1,1}(1,1,:)));
ci_stp = 1; % step size, multiply by 2 to get milliseconds
numboots = 500; % number of resamples

% loop through channels
for ci_ch_ct = 1:1:numchans
    for ci_tp = 1:ci_stp:tpmua2 % loop through time points
        MUAtmp = squeeze(MUA(ci_ch_ct,:,ci_tp)); % get the data from a time point
        ci_tmp = bootci(numboots,@mean,MUAtmp'); % get bootstrap derived CIs
        cimean(ci_ch_ct,ci_tp) = mean(MUAtmp); % record mean
        cilo(ci_ch_ct,ci_tp) = ci_tmp(1,1); % lower CI
        cihi(ci_ch_ct,ci_tp) = ci_tmp(2,1); % upper CI
        
        if cilo(ci_ch_ct,ci_tp)>0
            sigtimes(ci_ch_ct,ci_tp) = 1;
        elseif cihi(ci_ch_ct,ci_tp)<0
            sigtimes(ci_ch_ct,ci_tp) = -1;
        else
            sigtimes(ci_ch_ct,ci_tp) = 0;
        end
    end
end

figure
plotstp = ci_stp*2;
for ciplotct = 1:1:numchans
    subplot(numchans,1,ciplotct)
    plot(-250:plotstp:tpmua2,cihi(ciplotct,1:ci_stp:tpmua2))
    hold on
    plot(-250:plotstp:tpmua2,cilo(ciplotct,1:ci_stp:tpmua2))
    hold on
    plot(-250:plotstp:tpmua2,cimean(ciplotct,1:ci_stp:tpmua2))
end

figure
plotstp = ci_stp*2;
increases=[0:0.7:numchans*0.7];
for ciplotct = 1:1:numchans
    
    plot(-250:plotstp:tpmua2,cihi(ciplotct,1:ci_stp:tpmua2)+increases(ciplotct),'r')
    hold on
    plot(-250:plotstp:tpmua2,cilo(ciplotct,1:ci_stp:tpmua2)+increases(ciplotct),'r')
    hold on
    plot(-250:plotstp:tpmua2,cimean(ciplotct,1:ci_stp:tpmua2)+increases(ciplotct),'b')
    hold on
    plot(-250:plotstp:tpmua2,(zeros(length(-250:plotstp:tpmua2))+increases(ciplotct)),'k')
    
end
ylim([-0.5,15])

%% Cross correlation analysis
% now we have multiple simultaneous recordings and want to see how they are
% functionally connected 
numchans = length(all{1,1}(:,1,1));
numtrs = length(all{1,1}(1,:,1));
numtps = length(squeeze(all{1,1}(1,1,:)));

% making a meshgrid for looping through all combos of channels
chs1=1:1:numchans;
chs2=1:1:numchans;
[chs1,chs2]=meshgrid(chs1,chs2);
chs1=chs1(:);
chs2=chs2(:);
chcombos = [chs1,chs2];



corrs = {zeros(numtrs,length(chs1))};
xlags = {zeros(numtrs,length(chs1))};
negpeaks = zeros(numtrs,length(chs1));
negpeaklags = zeros(numtrs,length(chs1));
pospeaks = zeros(numtrs,length(chs1));
pospeaklags = zeros(numtrs,length(chs1));

% cross correlate all channel combos for every trial
for chanct = 1:length(chs1)
    
    for tr_ct = 1:numtrs
        
        
        xtmp = squeeze(site1LFP(chs1(chanct),tr_ct,:)); % one trial from site 1 (e.g. MGB)
        ytmp = squeeze(site2CSD(chs2(chanct),tr_ct,:)); % same trial from site 2
        [r,lags] = xcorr(xtmp,ytmp,'coeff'); %cross correlation
        corrs{tr_ct,chanct} = r;
        xlags{tr_ct,chanct} = lags';
        
        
        
        % recording the peaks in the xcorrelogram
        if min(r)<0
            negpeaks(tr_ct,chanct) = min(r(r<0));
            negpeaklags(tr_ct,chanct) = lags(r==negpeaks(tr_ct,chanct));
        else
            negpeaks(tr_ct,chanct) = NaN;
            negpeaklags(tr_ct,chanct) = NaN;
        end
        
        if max(r)>0
            pospeaks(tr_ct,chanct) = max(r(r>0));
            pospeaklags(tr_ct,chanct) = lags(r==pospeaks(tr_ct,chanct));
        else
            pospeaks(tr_ct,chanct) = NaN;
            pospeaklags(tr_ct,chanct) = NaN;
        end
        
    end
    
end

%% cross correlate all channel combos for avg resp

corrs = {zeros(length(chs1))};
xlags = {zeros(length(chs1))};
negpeaks = zeros(length(chs1));
negpeaklags = zeros(length(chs1));
pospeaks = zeros(length(chs1));
pospeaklags = zeros(length(chs1));


for chanct = 1:length(chs1)

       
        xtmp = squeeze(mean(site1LFP(chs1(chanct),:,:),2)); % average resp across trials
        ytmp = squeeze(mean(site2CSD(chs2(chanct),:,:),2)); % same ch combo from second site
        [r,lags] = xcorr(xtmp,ytmp,'coeff'); %cross correlation
        corrs{chanct,1} = r;
        xlags{chanct,1} = lags';
        
        
        
        % recording the peaks in the xcorrelogram
        if min(r)<0
            negpeaks(chanct) = min(r(r<0));
            negpeaklags(chanct) = lags(r==negpeaks(chanct));
        else
            negpeaks(chanct) = NaN;
            negpeaklags(chanct) = NaN;
        end
        
        if max(r)>0
            pospeaks(chanct) = max(r(r>0));
            pospeaklags(chanct) = lags(r==pospeaks(chanct));
        else
            pospeaks(chanct) = NaN;
            pospeaklags(chanct) = NaN;
        end
        
end
   


%% select channels of interest that get sorted/plotted
chcombo1 = 146;
chcombo2 = 348;
chcombo3 = 187;
chcombo4 = 334;
% chcombo1 = 70;
% chcombo2 = 327;
% chcombo3 = 84;
% chcombo4 = 147;

alt1 = 70;
alt3 = 84;
alt4 = 147;
%%
% plot a few trials' cross correlations
figure
legend1=('MGD -> supra. A1');
legend2=('MGV -> gran. A1');
for sorttrs = 20:1:22
    
    subplot(2,3,sorttrs-19)
    plot(xlags{sorttrs,chcombo1}(:),corrs{sorttrs,chcombo1}(:))
    hold on 
    plot(xlags{sorttrs,chcombo2}(:),corrs{sorttrs,chcombo2}(:))
    ylim([-0.6,0.6])
    ylabel('Cross Corr. Coeff')
    xlabel(' lag (ms)')
    title(['MGB -> A1 cross corr. Trial #', sprintf(num2str(sorttrs))])
    legend(legend1,legend2)
    set(gca,'fontsize', 16) 
    
end

subplot(2,3,4)
histogram(negpeaklags(:,chcombo1),'BinWidth',10,'Normalization','pdf') 
hold on
histogram(negpeaklags(:,chcombo2),'BinWidth',10,'Normalization','pdf')
xlabel('X corr peak lag')
ylabel('Proportion')
legend(legend1,legend2)
title('Negative X corr peaks')
set(gca,'fontsize', 16) 

hold on
subplot(2,3,6)
histogram(pospeaklags(:,chcombo1),'BinWidth',10,'Normalization','pdf') 
set(gca,'fontsize', 16) 
hold on
histogram(pospeaklags(:,chcombo2),'BinWidth',10,'Normalization','pdf') 
xlabel('X corr peak lag')
ylabel('Proportion')
legend(legend1,legend2)
title('Positive X corr peaks')

set(gca,'fontsize', 16) 

%% plotting x corr of avg resp
figure
subplot(1,2,1)
plot(xlags{chcombo1}(:),corrs{chcombo1}(:))
hold on
plot(xlags{chcombo2}(:),corrs{chcombo2}(:))
ylim([-0.8,0.8])
ylabel('Cross Corr. Coeff')
xlabel('lag(ms)')
title(['MGB -> A1 cross corr.'])
set(gca,'fontsize', 16)
legend('MGD -> supra. A1', 'MGV -> gran. A1')

subplot(1,2,2)
plot(xlags{chcombo3}(:),corrs{chcombo3}(:))
hold on
plot(xlags{chcombo4}(:),corrs{chcombo4}(:))
ylim([-0.8,0.8])
ylabel('Cross Corr. Coeff')
xlabel(' lag (ms)')
title(['MGB -> A1 cross corr.'])
set(gca,'fontsize', 16)
legend('MGD -> infra. A1', 'MGV -> infra. A1')


%% stats
xpos=pospeaklags(:,chcombo1);
ypos=pospeaklags(:,chcombo2);
xneg=negpeaklags(:,chcombo1);
yneg=negpeaklags(:,chcombo2);

% test for difference in distributions
[k_neg,p_neg] = kstest2(xneg,yneg)
[k_pos,p_pos] = kstest2(xpos,ypos)
