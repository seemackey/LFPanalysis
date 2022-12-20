%% Loops through the ephys extract function and gets out LFP data

clear
close all

%% path info 


% paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu009010\1-bu009010034@os';
%    '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\contproc\bu009010\2-bu009010034@os'};

paths = {'C:\Users\cmackey\Documents\MATLAB\1-bu009010034@os';
    'C:\Users\cmackey\Documents\MATLAB\2-bu009010034@os'
    };

%% loop through the function to extract data and find outliers
all=cell(length(paths));

for loopct = 1:1:length(paths)
    
    [CSD, LFP, trig] =  EphysExtractFxn(paths{loopct,1}); % extract data and stim triggers

    %unclean data
    all{1,loopct} = num2cell(CSD); 
    all{2,loopct} = num2cell(LFP);
    all{3,loopct} = num2cell(trig);
    
    [~, ~, ~, outlieridxtmp] = rejectartifacts(CSD, LFP, trig); % finding artifacts
    
    outliers{:,loopct} = outlieridxtmp; %at the end this has idx of outliers from both recordings

    
end
         

%% remove outliers across two sites (only works for two right now)
outliers = [outliers{:,1};outliers{:,2}];
outliersunique = unique(outliers);

for outlierct = 1:length(paths)
  
    all{1,outlierct}(:,outliersunique,:) = []; % CSD
    all{2,outlierct}(:,outliersunique,:) = []; % LFP
    all{3,outlierct}(:,outliersunique,:) = []; % triggers
    
end
    

%% Cross correlation analysis
% now we have multiple simultaneous recordings and want to see how they are
% functionally connected 
plott = 0;
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

site1LFP = cell2mat(all{2,1}(:,:,:)); %using LFP for site 1 (mgb)
site2CSD = cell2mat(all{1,2}(:,:,:)); % 2nd site, ctx

corrs = {zeros(numtrs,length(chs1))};
xlags = {zeros(numtrs,length(chs1))};

% cross correlate all channel combos for every trial
for chanct = 1:length(chs1)
    
    for tr_ct = 1:numtrs

        xtmp = squeeze(site1LFP(chs1(chanct),tr_ct,:)); % one trial from site 1 (e.g. MGB)
        ytmp = squeeze(site2CSD(chs2(chanct),tr_ct,:)); % same trial from site 2
        [r,lags] = xcorr(xtmp,ytmp,'coeff');
        corrs{tr_ct,chanct} = r;
        xlags{tr_ct,chanct} = lags';
        
    end
    
end


%% sorting


% check CSDs if needed
if plott == 1
    figure
    axpos=[0.1 0.1 0.8 0.8];
    figureax1a=axes('Position',axpos);
    [cax2] = csd_maker_no_subplot07(squeeze(mean(site2CSD(:,:,:),2)),(-200:2:200),1,[-10 200],[0 0],[],axpos,figureax1a);
    colormap(figureax1a,flipud(jet))
end

% plot a few trials' cross correlations
figure
for sorttrs = 30:1:32
    subplot(1,3,sorttrs-29)
    plot(xlags{sorttrs,45}(:),corrs{sorttrs,45}(:))
    hold on 
    plot(xlags{sorttrs,55}(:),corrs{sorttrs,55}(:))
    ylim([-0.5,0.5])
    ylabel('Cross Corr. Coeff')
    xlabel('Cross Corr. lag Re: noise onset(ms)')
    title(['MGB -> A1 cross corr. Trial #', sprintf(num2str(sorttrs))])
    legend('MGD -> A1 layer 1','MGV-> A1 layer 4')
end

