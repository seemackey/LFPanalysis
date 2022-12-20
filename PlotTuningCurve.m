clear all
close all

%% path info 

dir1 = 'C:\Users\cmackey\Documents\MATLAB';
paths = {'\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\bu015016\bu015016_tono\1-bu015016029@ep1';
   '\\NKI-LAKATOSLAB\lakatoslab_alpha\Buster\bu015016\bu015016_tono\2-bu015016029@ep1'};
fname='bu015016029_MGBmd';

load(paths{1,1}); % load the first file

%% some inputs, don't change
newadrate = 1000;
epoch_tframe            = [-50 100];
filtere = [0.5 300];
filteru = [300 5000];
baseline    = [-15 0];
seltime = [10 40];
disptime    = [-50 100];
 filtertype=1;
colors1=[0	0	0.5
    0	0	1
    0	0.5	1
    0	0.5	0.5
    0	0.5	0
    0	1	0
    0.5	1	0
    0.9	0.8	0
    1	0.5	0
    0.5	0.5	0
    0.5	0	0
    0.5	0	0.5
    0.5	0	1
    1	0	1
    1	0	0];
linecolor1=[0.3 0.3 0.3];
freqxaxis = [0.35,0.5,0.7,1,1.4,2,2.8,4,5.6,8,11.3,16,22.6,32];

%% epoching and filtering
time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);
[cnt_arej, cnte, cntm, cntc, cntu, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);
[trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig, params);

trig0=trig1s{1};
ttype1=ttype1s{1};% the order of the diff trig codes representing different freq tones
 
 
%% trigger downsample
    for i1=1:length(trig0)
        trig01(i1)    = round(trig0(i1)./(trig.adrate/newadrate)); % times of all trigs
    end
%      %filter the cnt in 1-100Hz to make nice CSD picts
%                      %to bandpass filter in the 1.2-2.5Hz range the MUA (cntm)
%                      n = 2;
%                      Wn = [1 100]/(newadrate/2); %the first 2 numbers are the freq range you want to pass, which is =/- 20% of the reprate
%                      [b,a] = butter(n,Wn);
%                      for i=1:size(cntc,1) %this is for the electrode chs
%                          cntc_ff(i,:)=filtfilt(b,a,cntc(i,:)); %creates a new variable which contains filtered data
%                          cntm_ff(i,:)=filtfilt(b,a,cntm(i,:));
%                      end
                 
    %%% to make epochs of the CSD and MUA
        x1 = round(epoch_tframe(1)*(newadrate/1000));
        x2 = round(epoch_tframe(2)*(newadrate/1000));
        for i=1:length(trig01)
            eegc(:,i,:)=cntc(:,trig01(i)+x1:trig01(i)+x2); %epoching the csd
        end
        for i=1:length(trig01)
            eegm(:,i,:)=cntm(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
        end
%% Artifact reject
            cmax = [];
            for i1=1:size(eegc,2)
                cmax(i1)=max(mean(abs(squeeze(eegc(:,i1,:))),1));
            end
            [b,idxC,outliers] = deleteoutliers(cmax,0.5);
            mmax = [];
            for i1=1:size(eegm,2)
                mmax(i1)=max(mean(abs(squeeze(eegm(:,i1,:))),1));
            end
            [b,idxM,outliers] = deleteoutliers(mmax,0.5);
            idx=[idxC idxM];
            idx = unique(idx);
            trig01(:,idx)=[];
            ttype1(idx,:)=[];
            eegc(:,idx,:)           = []; % gets rid of the epochs with artifacts
            eegm(:,idx,:)           = []; % gets rid of the epochs with artifacts
            
            
%%            %triggers for tono
            a=0; trigtype=[];
            for i=1:256
                z=find(ttype1==i);
                if ~isempty(z)
                    a=a+1;
                    trigtype(a)=i;
                end
            end
            %sweepnos=size(eegm(:,ttype1==trigtype(1),:),2);
            
              %baseline
          for i=1:size(eegm,1)
            for iii=1:size(eegm,2)
                eegm(i,iii,:)=squeeze(eegm(i,iii,:))-squeeze(mean(eegm(i,iii,max(find(time<=baseline(1))):max(find(time<=baseline(2)))),3));
            end
          end
  
        XX= squeeze(mean(mean(eegm(13:15,:,max(find(time<=seltime(1))):max(find(time<=seltime(2)))),1),3)); %chs 13=15 are gran layer chs for file 1-gt005013
%         XX1= squeeze(mean(eegm(:,:,max(find(time<=seltime(1))):max(find(time<=seltime(2)))),3)); %avgs across all chs
         % smoothing after measuring amplitude, for display
        for i=1:size(eegm,1)
            for iii=1:size(eegm,2)
                for times_smooth=1:1
                    eegm(i,iii,:)=smooth(squeeze(eegm(i,iii,:)),'moving',5);
                end
            end
        end
        
            %For avglam layer
            
                for trigcik=1:length(trigtype)
                    z = find(ttype1==trigtype(trigcik) );
                    amps{trigcik}           = XX(z);
                    tunamp_mean(trigcik)    = mean(XX(z),2);
                    tunamp_stderr(trigcik)  = std(XX(z),0,2)/sqrt(length(z));
                    avgs(trigcik,:)=squeeze(mean(mean(eegm(:,z,:),1),2));
                end
               
                
                
                y=tunamp_mean(1,1:14); %this is the amp value of each frq
                x=[353 500 706 1000 1414 2000 2828 4000 5656 8000 11312 16000 22614 32000];
                y1gran=max(y);
                y = y /y1gran;%normalize to max val
                N = length(y);
                lev50 = 0.5;
                [r, BF] = find(y==max(y)); %BF is best freq of the gran layer
                tunamp_mean1=tunamp_mean./y1gran;%normaisining by max MUA value of BF
                BFgran=x(BF);
                 clear z
                 
                 z = find(ttype1==trigtype(BF) );
%% plotting CSD
                 
                 curfig = figure;
                 set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                 axpos=[0.1 0.1 0.8 0.8];
                 figureax1a=axes('Position',axpos);
                 [cax2] = csd_maker_no_subplot07(squeeze(mean(eegc(:,z,:),2)),time,1,[-10 100],[0 0],[],axpos,figureax1a);
                 colormap(figureax1a,flipud(jet))
                 title (['BF= ' num2str(x(BF)) ',time interval: ' num2str(seltime(1)) ' - ' num2str(seltime(2)),])
                 axes('Position',[0 0.98 1 0.2],'Visible','off');
                 text(0.5,0,[fname '  BF CSD laminar response, scale: ' num2str(cax2(1),'%11.3g') '-' num2str(cax2(2),'%11.3g')],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                 print('-dpsc','-r1200','-bestfit',[dir1  fname '  BF CSD laminar response' '.ps'])
                 print('-djpeg','-r300',[dir1 fname ' BF CSD laminar response' '.jpg'])
                 
                 curfig = figure;
                     set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                subplot(2,1,1)
                plot(1:size(tunamp_mean(1,1:14),2),tunamp_mean1(1,1:14),'color','b')
                set(gca,'xlim',[0.5 size(tunamp_mean(1,1:14),2)+0.5])
                set(gca,'xtick',1:1:size(tunamp_mean(1,1:14),2),'xticklabel',trigtype)
                set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',8)
                title (['Tuning curve' ])
                 subplot(2,1,2)
                for i1=1:size(avgs,1)
                    plot(time,avgs(i1,:),'color',colors1(i1,:),'linewidth',1)
                    hold on
                end
                 title (['granMUA of the 14 tones & BBN overlaid' ])
                  axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,[fname ' Tuning curve & MUA resp to tone' ],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                     print('-dpsc','-r1200','-bestfit',[dir1  fname '   Tuning curve & MUA resp to tone' '.ps'])
                    print('-djpeg','-r300',[dir1 fname '  Tuning curve & MUA resp to tone' '.jpg'])
                
                 tunamp_meanGT005=tunamp_mean1;
                    avgsGT005=avgs;
                    
                    
% plot MUA   
%smoothing for display
for i=1:size(eegm,1)
 for iii=1:size(eegm,2)
     for times_smooth=1:2
         eegmbS(i,iii,:)=smooth(squeeze(eegm(i,iii,:)),'moving',5);%smoothed data
     end
 end
end


curfig = figure;
set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')

axpos=[0.1 0.1 0.8 0.8];
%                     axpos=[0.1 0.4 0.5 0.5];
figureax1a=axes('Position',axpos);
[cax2] = csd_maker_no_subplot07(squeeze(mean(eegm(:,:,:),2)),time,0,[-10 100],[0 0],[0 0],axpos,figureax1a);%MUA make 0 
colormap(figureax1a,hot)
title([ fname 'Tone Evoked MUA,caxis ' num2str(cax2)  ',N=' num2str(length(trig01)) ],'FontSize',9)



%% clear variables
clear   avgs eegm eegc cntm cntc BF BFgran amps trig trig01 trig0 ttype1 tunampmean XX



%% load A1 file for comparison %% 
load(paths{2,1}) %load second file
fname='bu015016029_A1';
%

                    [cnt_arej, cnte, cntm, cntc, cntu, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);
                    [trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig, params);
                    
                    trig0=trig1s{1};
                    ttype1=ttype1s{1};% the order of the diff trig codes representing different freq tones
                    %%%%%%%%% trigger downsample
                    for i1=1:length(trig0)
                        trig01(i1)    = round(trig0(i1)./(trig.adrate/newadrate)); % times of all trigs
                    end
                    %      %filter the cnt in 1-100Hz to make nice CSD picts
                    %                      %to bandpass filter in the 1.2-2.5Hz range the MUA (cntm)
                    %                      n = 2;
                    %                      Wn = [1 100]/(newadrate/2); %the first 2 numbers are the freq range you want to pass, which is =/- 20% of the reprate
                    %                      [b,a] = butter(n,Wn);
                    %                      for i=1:size(cntc,1) %this is for the electrode chs
                    %                          cntc_ff(i,:)=filtfilt(b,a,cntc(i,:)); %creates a new variable which contains filtered data
                    %                          cntm_ff(i,:)=filtfilt(b,a,cntm(i,:));
                    %                      end
                    
                    %%% to make epochs of the CSD and MUA
                    x1 = round(epoch_tframe(1)*(newadrate/1000));
                    x2 = round(epoch_tframe(2)*(newadrate/1000));
                    for i=1:length(trig01)
                        eegc(:,i,:)=cntc(:,trig01(i)+x1:trig01(i)+x2); %epoching the csd
                    end
                    for i=1:length(trig01)
                        eegm(:,i,:)=cntm(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
                    end
                    %Artifact reject
                    cmax = [];
                    for i1=1:size(eegc,2)
                        cmax(i1)=max(mean(abs(squeeze(eegc(:,i1,:))),1));
                    end
                    [b,idxC,outliers] = deleteoutliers(cmax,0.5);
                    mmax = [];
                    for i1=1:size(eegm,2)
                        mmax(i1)=max(mean(abs(squeeze(eegm(:,i1,:))),1));
                    end
                    [b,idxM,outliers] = deleteoutliers(mmax,0.5);
                    idx=[idxC idxM];
                    idx = unique(idx);
                    trig01(:,idx)=[];
                    ttype1(idx,:)=[];
                    eegc(:,idx,:)           = []; % gets rid of the epochs with artifacts
                    eegm(:,idx,:)           = []; % gets rid of the epochs with artifacts
                    
%% extract tone resps for tuning curve 
                    a=0; trigtype=[];
                    for i=1:256
                        z=find(ttype1==i);
                        if ~isempty(z)
                            a=a+1;
                            trigtype(a)=i;
                        end
                    end
                    %sweepnos=size(eegm(:,ttype1==trigtype(1),:),2);
                    
                    %baseline
                    for i=1:size(eegm,1)
                        for iii=1:size(eegm,2)
                            eegm(i,iii,:)=squeeze(eegm(i,iii,:))-squeeze(mean(eegm(i,iii,max(find(time<=baseline(1))):max(find(time<=baseline(2)))),3));
                        end
                    end
                    
                    XX= squeeze(mean(mean(eegm(12:13,:,max(find(time<=seltime(1))):max(find(time<=seltime(2)))),1),3)); %chs 12-13 are gran layer chs for file 1-gt006015
                    %         XX1= squeeze(mean(eegm(:,:,max(find(time<=seltime(1))):max(find(time<=seltime(2)))),3)); %avgs across all chs
                    % smoothing after measuring amplitude, for display
                    for i=1:size(eegm,1)
                        for iii=1:size(eegm,2)
                            for times_smooth=1:1
                                eegm(i,iii,:)=smooth(squeeze(eegm(i,iii,:)),'moving',5);
                            end
                        end
                    end
                    
                    %For avglam layer
                    for trigcik=1:length(trigtype)
                        z = find(ttype1==trigtype(trigcik) );
                        amps{trigcik}           = XX(z);
                        tunamp_mean(trigcik)    = mean(XX(z),2);
                        tunamp_stderr(trigcik)  = std(XX(z),0,2)/sqrt(length(z));
                        avgs(trigcik,:)=squeeze(mean(mean(eegm(:,z,:),1),2));
                    end
                    
                   
                    y=tunamp_mean(1,1:14); %this is the amp value of each frq
                    x=[353 500 706 1000 1414 2000 2828 4000 5656 8000 11312 16000 22614 32000];
                    y1gran=max(y);
                    y = y /y1gran;%normalize to max val
                    N = length(y);
                    lev50 = 0.5;
                    [r, BF] = find(y==max(y)); %BF is best freq of the gran layer
                     tunamp_mean1=tunamp_mean./y1gran; %normaisining by max MUA value of BF
                    BFgran=x(BF);
                    clear z
                    
                    z = find(ttype1==trigtype(BF) );
                    
                    curfig = figure;
                    set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                    axpos=[0.1 0.1 0.8 0.8];
                    figureax1a=axes('Position',axpos);
                    [cax2] = csd_maker_no_subplot07(squeeze(mean(eegc(:,z,:),2)),time,1,[-10 100],[0 0],[],axpos,figureax1a);
                    colormap(figureax1a,flipud(jet))
                    title (['BF= ' num2str(x(BF)) ',time interval: ' num2str(seltime(1)) ' - ' num2str(seltime(2)),])
                    axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,[fname '  BF CSD laminar response, scale: ' num2str(cax2(1),'%11.3g') '-' num2str(cax2(2),'%11.3g')],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                    print('-dpsc','-r1200','-bestfit',[dir1  fname '  BF CSD laminar response' '.ps'])
                    print('-djpeg','-r300',[dir1 fname ' BF CSD laminar response' '.jpg'])
                    
                    curfig = figure;
                    set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                    subplot(2,1,1)
                    plot(1:size(tunamp_mean(1,1:14),2),tunamp_mean(1,1:14),'color','b')
                    set(gca,'xlim',[0.5 size(tunamp_mean(1,1:14),2)+0.5])
                    set(gca,'xtick',1:1:size(tunamp_mean(1,1:14),2),'xticklabel',trigtype)
                    set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',8)
                    title (['Tuning curve' ])
                    subplot(2,1,2)
                    for i1=1:size(avgs,1)
                        plot(time,avgs(i1,:),'color',colors1(i1,:),'linewidth',1)
                        hold on
                    end
                    title (['granMUA of the 14 tones & BBN overlaid' ])
                    axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,[fname ' Tuning curve & MUA resp to tone' ],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                    print('-dpsc','-r1200','-bestfit',[dir1  fname '   Tuning curve & MUA resp to tone' '.ps'])
                    print('-djpeg','-r300',[dir1 fname '  Tuning curve & MUA resp to tone' '.jpg'])
                    
                    tunamp_meanGT006=tunamp_mean1;
                    avgsGT006=avgs;
                    
                    %normalising the avg MUA by max MUA value
                    Xgt005=avgsGT005(13,:)./max(avgsGT005(13,:));
                    Xgt006=avgsGT006(13,:)./max(avgsGT006(13,:));
                    
                   
                   

                    % plot tuning curves - Chase
                    
                    curfig = figure;
                    set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                    subplot(2,1,1)
                    semilogx(freqxaxis,tunamp_meanGT005(1,1:14),'color','b')
                    hold on
                    semilogx(freqxaxis,tunamp_meanGT006(1,1:14),'color','r')
                    %set(gca,'xtick',1:1:size(tunamp_mean(1,1:15),2),'xticklabel',trigtype)
                     subplot(2,1,2)
                     plot(time,Xgt005,'color','b','linewidth',1)
                     hold on
                      plot(time,Xgt006,'color','r','linewidth',1)
                     axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,[ 'Comparison of Gt006(A1,red) & Gt005(belt,blue)tuning, BF=22KHz' ],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                    print('-dpsc','-r1200','-bestfit',[dir1  'Comparison of Gt006(A1,red) & Gt005(belt,blue)tuning, BF=22KHz' '.ps'])
                    print('-djpeg','-r300',[dir1 'Comparison of Gt006(A1,red) & Gt005(belt,blue)tuning, BF=22KHz' '.jpg'])
                    
                    
                                        %% plot MUA   
%smoothing for display
for i=1:size(eegm,1)
 for iii=1:size(eegm,2)
     for times_smooth=1:2
         eegmbS(i,iii,:)=smooth(squeeze(eegm(i,iii,:)),'moving',5);%smoothed data
     end
 end
end


curfig = figure;
set(curfig,'position',[100   50  800   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')

axpos=[0.1 0.1 0.8 0.8];
%                     axpos=[0.1 0.4 0.5 0.5];
figureax1a=axes('Position',axpos);
[cax2] = csd_maker_no_subplot07(squeeze(mean(eegm(:,:,:),2)),time,0,[-10 100],[0 0],[0 0],axpos,figureax1a);%MUA make 0
colormap(figureax1a,hot)
title([ fname 'Tone Evoked MUA,caxis ' num2str(cax2)  ',N=' num2str(length(trig01)) ],'FontSize',9)