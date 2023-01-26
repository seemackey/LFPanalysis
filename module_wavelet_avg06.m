function [wraw,wlfp,wai,epoch_tframe] = module_wavelet_avg06(wraw, wlfp, wai, trig, params, epoch_tframe, bool)

% epoch_tframe = [-50 150];
epoch_tframe = epoch_tframe;

[trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig, params);

wraw.avg.poe          = {};
wraw.avg.pom          = {};
wraw.avg.pob          = {};
wraw.avg.poc          = {};
wraw.avg.poeye        = {};
wraw.avg.posurface  = {};

wraw.avg.itce          = {};
wraw.avg.itcm          = {};
wraw.avg.itcb          = {};
wraw.avg.itcc          = {};
wraw.avg.itceye        = {};
wraw.avg.itcsurface  = {};

wraw.avg.ame          = {};
wraw.avg.amm          = {};
wraw.avg.amb          = {};
wraw.avg.amc          = {};
wraw.avg.ameye        = {};
wraw.avg.amsurface  = {};

trigcik_counter = 0;

adrate  = wraw.adrate;
dt      = 1000/adrate;

for trigcik = 1:length(trig1s)
    
    trig1           = trig1s{trigcik};
    ttype1          = ttype1s{trigcik};
    triglength1     = triglength1s{trigcik};
    findex2         = findex2s{trigcik};
    
    time            = dt*round(epoch_tframe(1)/dt):dt:dt*round(epoch_tframe(2)/dt);
    
    wraw.avg.time = time;
    
    if size(ttype1,2)>size(ttype1,1)
        ttype1=ttype1';
    end
    
    if ~isempty(wraw.adrate)
        x1 = ceil(epoch_tframe(1)/dt);
        x2 = floor(epoch_tframe(2)/dt);
        z=find(trig1<abs(x1) | trig1>length(wraw.cntc)-x2);
        trig1(z)        = [];
        
        ttype1(z)       = [];
        triglength1(z)  = [];
        
    elseif ~isempty(wlfp.adrate)
        x1 = ceil(epoch_tframe(1)/dt);
        x2 = floor(epoch_tframe(2)/dt);
        z=find(trig1<abs(x1) | trig1>length(wlfp.cnt)-x2);
        trig1(z)        = [];
        
        ttype1(z)       = [];
        triglength1(z)  = [];
        
    elseif ~isempty(wai.adrate)
        x1 = ceil(epoch_tframe(1)/dt);
        x2 = floor(epoch_tframe(2)/dt);
        z=find(trig1<abs(x1) | trig1>length(wai.cnt)-x2);
        trig1(z)        = [];
        
        ttype1(z)       = [];
        triglength1(z)  = [];
    end
    
    if bool.arej == 1
        
        x1 = round(epoch_tframe(1)/dt);
        x2 = round(epoch_tframe(2)/dt);
        
        idx1=[];
        idx2=[];
        if ~isempty(wraw.cntc)
            for i1=1:length(trig1)
                eeg_po(:,i1,:,:)           = wraw.cntc_po(:,:,trig1(i1)+x1:trig1(i1)+x2);
            end
            mmax = [];
            for i1=1:size(eeg_po,2)
                mmax(i1)=max(squeeze(mean(mean(eeg_po(:,i1,:,:),1),3)));
            end
            [b,idx1,outliers] = deleteoutliers(mmax,0.05);
        end
        
        if ~isempty(wraw.cntm)
            for i1=1:length(trig1)
                eeg_po(:,i1,:,:)           = wraw.cntm_po(:,:,trig1(i1)+x1:trig1(i1)+x2);
            end
            mmax = [];
            for i1=1:size(eeg_po,2)
                mmax(i1)=max(squeeze(mean(mean(eeg_po(:,i1,:,:),1),3)));
            end
            [b,idx2,outliers] = deleteoutliers(mmax,0.05);
        end
        
        idx=[idx1 idx2];
        idx = unique(idx);
    else
        idx = [];
    end
    
    trig1(idx)        = [];
    ttype1(idx)       = [];
    triglength1(idx)  = [];
    
    wraw.avg.rejsweeps = idx;
    wlfp.avg.rejsweeps = idx;
    wai.avg.rejsweeps = idx;
    
    %%%% for pattern detection
    %     ttd         = diff(ttype1);
    %     ttd         = [0; ttd];
    %     z           = find(ttd>400 & ttype1>999 & ttype1<1999);
    %     ttype1(z)   = 666;
    
    z=find(ttype1>999 & ttype1<1999);
    ttype1(z)   = ttype1(z)-1000;
    
    if length(findex2)>3
        if  strcmp(findex2(1:3),'_Sp') && trigcik >1
            ttype1(1:end) = 666;
        end
    end
    
    if ~isempty (trig1)
        
        trigcik_counter = trigcik_counter+1;
        
        if ~isempty(wraw.adrate)
            
            
            a=0; trigtype=[];
            for i=0:20000
                z=find(ttype1==i);
                if ~isempty(z)
                    a=a+1;
                    trigtype(a)=i;
                end
            end
            
            avgtrigtype    = [];
            avgsweepno     = [];
            
            if length(trigtype)>1 && length(trigtype)<4
                % if trigtype(2)==2 || trigtype(3)==12
                    
                    z=find(ttype1==2 | ttype1==12);
                    ttype1(z)=10001;
                    
                    z=[0;z; length(ttype1)];
                    for i1=1:length(z)-1
                        a=0;
                        for i2=z(i1)+1:z(i1+1)-1
                            a=a+1;
                            ttype1(i2)=a;
                        end
                    end
                    if ttype1(end)<10001
                        ttype1(end)=ttype1(end-1)+1;
                    end
                    
                    %%%% have to redo trigger types
                    a=0; trigtype=[];
                    for i=0:20000
                        z=find(ttype1==i);
                        if ~isempty(z)
                            a=a+1;
                            trigtype(a)=i;
                        end
                    end
                % end
            end
            
            
            
            if ~isempty(wraw.cntc)
                
                for trigcik2 = 1:length(trigtype)+1;
                    
                    if trigcik2 == 1
                        z = 1:length(ttype1);
                    else
                        z=find(ttype1==trigtype(trigcik2-1));
                    end
                    trig001 = trig1(z);
                    ttype001 = ttype1(z);
                    
                    x1 = round(epoch_tframe(1)/dt);
                    x2 = round(epoch_tframe(2)/dt);
                    
                    avgsweepno(trigcik2) = length(trig001);
                    if trigcik2 == 1
                        avgtrigtype(trigcik2) = -1;
                    else
                        avgtrigtype(trigcik2) = trigtype(trigcik2-1);
                    end
                    
                    eeg_po=zeros(size(wraw.cntc_po,1),length(trig001),size(wraw.cntc_po,2),length(time));
                    eeg_ph=zeros(size(wraw.cntc_po,1),length(trig001),size(wraw.cntc_po,2),length(time));
                    for i1=1:length(trig001)
                        eeg_po(:,i1,:,:)           = wraw.cntc_po(:,:,trig001(i1)+x1:trig001(i1)+x2);
                        eeg_ph(:,i1,:,:)           = wraw.cntc_ph(:,:,trig001(i1)+x1:trig001(i1)+x2);
                    end
                    
                    
                    
                    
                    if trigcik2 == 1
                        am = zeros(length(trigtype)+1,size(eeg_ph,1),size(eeg_ph,3),size(eeg_ph,4));
                        r = zeros(size(am));
                        p = zeros(size(am));
                        avgpo = zeros(size(am));
                    end
                    avgpo(trigcik2,:,:,:)=squeeze(mean(eeg_po(:,:,:,:),2));
                    for chancik=1:size(eeg_ph,1)
                        for i1=1:size(eeg_ph,3)
                            [p(trigcik2,chancik,i1,:), r(trigcik2,chancik,i1,:)] = rayleigh(squeeze(eeg_ph(chancik,:,i1,:)));
                            for i2=1:size(eeg_ph,4)
                                [am(trigcik2,chancik,i1,i2), ad] = angstat(squeeze(eeg_ph(chancik,:,i1,i2))');
                            end
                        end
                    end
                    trig001s{trigcik2} =  trig001;
                    ttype001s{trigcik2} =  ttype001;
                    
                end
                wraw.avg.sweepno{trigcik_counter}    = avgsweepno;
                wraw.avg.strigtype{trigcik_counter}  = avgtrigtype;
                wraw.avg.poc{trigcik_counter}        = avgpo;
                wraw.avg.itcc{trigcik_counter}        = r;
                wraw.avg.itcpc{trigcik_counter}        = p;
                wraw.avg.amc{trigcik_counter}        = am;
                wraw.avg.ttype{trigcik_counter}      = ttype001s;
                wraw.avg.trig{trigcik_counter}      = trig001s;
            end
            
            
            if ~isempty(wraw.cntm)
                
                for trigcik2 = 1:length(trigtype)+1;
                    
                    if trigcik2 == 1
                        z = 1:length(ttype1);
                    else
                        z=find(ttype1==trigtype(trigcik2-1));
                    end
                   trig001 = trig1(z);
                    ttype001 = ttype1(z);
                    
                    x1 = round(epoch_tframe(1)/dt);
                    x2 = round(epoch_tframe(2)/dt);
                    
                    avgsweepno(trigcik2) = length(trig001);
                    if trigcik2 == 1
                        avgtrigtype(trigcik2) = -1;
                    else
                        avgtrigtype(trigcik2) = trigtype(trigcik2-1);
                    end
                    
                    eeg_po=zeros(size(wraw.cntm_po,1),length(trig001),size(wraw.cntm_po,2),length(time));
                    eeg_ph=zeros(size(wraw.cntm_po,1),length(trig001),size(wraw.cntm_po,2),length(time));
                    for i1=1:length(trig001)
                        eeg_po(:,i1,:,:)           = wraw.cntm_po(:,:,trig001(i1)+x1:trig001(i1)+x2);
                        eeg_ph(:,i1,:,:)           = wraw.cntm_ph(:,:,trig001(i1)+x1:trig001(i1)+x2);
                    end
                    
                    
                    
                    if trigcik2 == 1
                        am = zeros(length(trigtype)+1,size(eeg_ph,1),size(eeg_ph,3),size(eeg_ph,4));
                        r = zeros(size(am));
                        p = zeros(size(am));
                        avgpo = zeros(size(am));
                    end
                    avgpo(trigcik2,:,:,:)=squeeze(mean(eeg_po(:,:,:,:),2));
                    for chancik=1:size(eeg_ph,1)
                        for i1=1:size(eeg_ph,3)
                            [p(trigcik2,chancik,i1,:), r(trigcik2,chancik,i1,:)] = rayleigh(squeeze(eeg_ph(chancik,:,i1,:)));
                            for i2=1:size(eeg_ph,4)
                                [am(trigcik2,chancik,i1,i2), ad] = angstat(squeeze(eeg_ph(chancik,:,i1,i2))');
                            end
                        end
                    end
                    
                     trig001s{trigcik2} =  trig001;
                    ttype001s{trigcik2} =  ttype001;
                    
                end
                
                wraw.avg.sweepno{trigcik_counter}    = avgsweepno;
                wraw.avg.strigtype{trigcik_counter}  = avgtrigtype;
                wraw.avg.pom{trigcik_counter}       = avgpo;
                wraw.avg.itcm{trigcik_counter}        = r;
                wraw.avg.itcpm{trigcik_counter}        = p;
                wraw.avg.amm{trigcik_counter}        = am;
                wraw.avg.ttype{trigcik_counter}      = ttype001s;
                wraw.avg.trig{trigcik_counter}      = trig001s;
                
            end
            
            if ~isempty(wraw.cntb)
                for trigcik2 = 1:length(trigtype)+1;
                    
                    if trigcik2 == 1
                        z = 1:length(ttype1);
                    else
                        z=find(ttype1==trigtype(trigcik2-1));
                    end
                    trig001 = trig1(z);
                    ttype001 = ttype1(z);
                    
                    x1 = round(epoch_tframe(1)/dt);
                    x2 = round(epoch_tframe(2)/dt);
                    
                    avgsweepno(trigcik2) = length(trig001);
                    if trigcik2 == 1
                        avgtrigtype(trigcik2) = -1;
                    else
                        avgtrigtype(trigcik2) = trigtype(trigcik2-1);
                    end
                    
                    eeg_po=zeros(size(wraw.cntb_po,1),length(trig001),size(wraw.cntb_po,2),length(time));
                    eeg_ph=zeros(size(wraw.cntb_po,1),length(trig001),size(wraw.cntb_po,2),length(time));
                    for i1=1:length(trig001)
                        eeg_po(:,i1,:,:)           = wraw.cntb_po(:,:,trig001(i1)+x1:trig001(i1)+x2);
                        eeg_ph(:,i1,:,:)           = wraw.cntb_ph(:,:,trig001(i1)+x1:trig001(i1)+x2);
                    end
                    
                    
                    
                    if trigcik2 == 1
                        am = zeros(length(trigtype)+1,size(eeg_ph,1),size(eeg_ph,3),size(eeg_ph,4));
                        r = zeros(size(am));
                        p = zeros(size(am));
                        avgpo = zeros(size(am));
                    end
                    avgpo(trigcik2,:,:,:)=squeeze(mean(eeg_po(:,:,:,:),2));
                    for chancik=1:size(eeg_ph,1)
                        for i1=1:size(eeg_ph,3)
                            [p(trigcik2,chancik,i1,:), r(trigcik2,chancik,i1,:)] = rayleigh(squeeze(eeg_ph(chancik,:,i1,:)));
                            for i2=1:size(eeg_ph,4)
                                [am(trigcik2,chancik,i1,i2), ad] = angstat(squeeze(eeg_ph(chancik,:,i1,i2))');
                            end
                        end
                    end
                     trig001s{trigcik2} =  trig001;
                    ttype001s{trigcik2} =  ttype001;
                end
                wraw.avg.sweepno{trigcik_counter}    = avgsweepno;
                wraw.avg.strigtype{trigcik_counter}  = avgtrigtype;
                wraw.avg.pob{trigcik_counter}        = avgpo;
                wraw.avg.itcb{trigcik_counter}        = r;
                wraw.avg.itcpb{trigcik_counter}        = p;
                wraw.avg.amb{trigcik_counter}        = am;
                wraw.avg.ttype{trigcik_counter}      = ttype001s;
                wraw.avg.trig{trigcik_counter}      = trig001s;
            end
            
            if ~isempty(wlfp.cnt)
                for trigcik2 = 1:length(trigtype)+1;
                    
                    if trigcik2 == 1
                        z = 1:length(ttype1);
                    else
                        z=find(ttype1==trigtype(trigcik2-1));
                    end
                   trig001 = trig1(z);
                    ttype001 = ttype1(z);
                    
                    x1 = round(epoch_tframe(1)/dt);
                    x2 = round(epoch_tframe(2)/dt);
                    
                    avgsweepno(trigcik2) = length(trig001);
                    if trigcik2 == 1
                        avgtrigtype(trigcik2) = -1;
                    else
                        avgtrigtype(trigcik2) = trigtype(trigcik2-1);
                    end
                    
                    eeg_po=zeros(size(wlfp.cnt_po,1),length(trig001),size(wlfp.cnt_po,2),length(time));
                    eeg_ph=zeros(size(wlfp.cnt_po,1),length(trig001),size(wlfp.cnt_po,2),length(time));
                    for i1=1:length(trig001)
                        eeg_po(:,i1,:,:)           = wlfp.cnt_po(:,:,trig001(i1)+x1:trig001(i1)+x2);
                        eeg_ph(:,i1,:,:)           = wlfp.cnt_ph(:,:,trig001(i1)+x1:trig001(i1)+x2);
                    end
                    
                    
                    if trigcik2 == 1
                        am = zeros(length(trigtype)+1,size(eeg_ph,1),size(eeg_ph,3),size(eeg_ph,4));
                        r = zeros(size(am));
                        p = zeros(size(am));
                        avgpo = zeros(size(am));
                    end
                    avgpo(trigcik2,:,:,:)=squeeze(mean(eeg_po(:,:,:,:),2));
                    for chancik=1:size(eeg_ph,1)
                        for i1=1:size(eeg_ph,3)
                            [p(trigcik2,chancik,i1,:), r(trigcik2,chancik,i1,:)] = rayleigh(squeeze(eeg_ph(chancik,:,i1,:)));
                            for i2=1:size(eeg_ph,4)
                                [am(trigcik2,chancik,i1,i2), ad] = angstat(squeeze(eeg_ph(chancik,:,i1,i2))');
                            end
                        end
                    end
                     trig001s{trigcik2} =  trig001;
                    ttype001s{trigcik2} =  ttype001;
                    
                end
                
                wlfp.avg.sweepno{trigcik_counter}    = avgsweepno;
                wlfp.avg.strigtype{trigcik_counter}  = avgtrigtype;
                wlfp.avg.po{trigcik_counter}       = avgpo;
                wlfp.avg.itc{trigcik_counter}        = r;
                wlfp.avg.itcp{trigcik_counter}        = p;
                wlfp.avg.am{trigcik_counter}        = am;
                wlfp.avg.ttype{trigcik_counter}      = ttype001s;
                wlfp.avg.trig{trigcik_counter}      = trig001s;
            end
            if ~isempty(wai.cnt)
                
                for trigcik2 = 1:length(trigtype)+1;
                    
                    if trigcik2 == 1
                        z = 1:length(ttype1);
                    else
                        z=find(ttype1==trigtype(trigcik2-1));
                    end
                    trig001 = trig1(z);
                    ttype001 = ttype1(z);
                    
                    x1 = round(epoch_tframe(1)/dt);
                    x2 = round(epoch_tframe(2)/dt);
                    
                    avgsweepno(trigcik2) = length(trig001);
                    if trigcik2 == 1
                        avgtrigtype(trigcik2) = -1;
                    else
                        avgtrigtype(trigcik2) = trigtype(trigcik2-1);
                    end
                    
                    eeg_po=zeros(size(wai.cnt_po,1),length(trig001),size(wai.cnt_po,2),length(time));
                    eeg_ph=zeros(size(wai.cnt_po,1),length(trig001),size(wai.cnt_po,2),length(time));
                    for i1=1:length(trig001)
                        eeg_po(:,i1,:,:)           = wai.cnt_po(:,:,trig001(i1)+x1:trig001(i1)+x2);
                        eeg_ph(:,i1,:,:)           = wai.cnt_ph(:,:,trig001(i1)+x1:trig001(i1)+x2);
                    end
                    
                    
                    
                    if trigcik2 == 1
                        am = zeros(length(trigtype)+1,size(eeg_ph,1),size(eeg_ph,3),size(eeg_ph,4));
                        r = zeros(size(am));
                        p = zeros(size(am));
                        avgpo = zeros(size(am));
                    end
                    avgpo(trigcik2,:,:,:)=squeeze(mean(eeg_po(:,:,:,:),2));
                    
                    for chancik=1:size(eeg_ph,1)
                        for i1=1:size(eeg_ph,3)
                            [p(trigcik2,chancik,i1,:), r(trigcik2,chancik,i1,:)] = rayleigh(squeeze(eeg_ph(chancik,:,i1,:)));
                            for i2=1:size(eeg_ph,4)
                                [am(trigcik2,chancik,i1,i2), ad] = angstat(squeeze(eeg_ph(chancik,:,i1,i2))');
                            end
                        end
                    end
                     trig001s{trigcik2} =  trig001;
                    ttype001s{trigcik2} =  ttype001;
                end
                
                wai.avg.sweepno{trigcik_counter}    = avgsweepno;
                wai.avg.strigtype{trigcik_counter}  = avgtrigtype;
                wai.avg.po{trigcik_counter}       = avgpo;
                wai.avg.itc{trigcik_counter}        = r;
                wai.avg.itcp{trigcik_counter}        = p;
                wai.avg.am{trigcik_counter}        = am;
                wai.avg.ttype{trigcik_counter}      = ttype001s;
                wai.avg.trig{trigcik_counter}      = trig001s;
            end
        end
    end
end


% curfig = figure;
% set(curfig,'position',[100   50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
% 
% t3 = wraw.avg.time;
% frq = wraw.frq;
% ylabels     = num2str(frq','%11.3g');
% colormap jet
% ylabelres=15;
% fsize = 6;
% 
% subplot(2,2,1)
% mapvar=squeeze(mean(wraw.avg.poc{1},1));
% surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
% set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
% set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
% caxset1=caxis;
% caxset1(2)=max(max(mapvar))*0.9;
% caxis(caxset1)
% title (['poc, ' num2str(wraw.avg.sweepno{1}(1)) ' epochs, sig:' num2str(rayleigh_p(wraw.avg.sweepno{1}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
% 
% subplot(2,2,2)
% mapvar=squeeze(mean(wraw.avg.itcc{1},1));
% surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
% set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
% set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
% caxset1=caxis;
% caxset1(2)=max(max(mapvar))*0.9;
% caxis(caxset1)
% title (['itcc, ' num2str(wraw.avg.sweepno{1}(1)) ' epochs, sig:' num2str(rayleigh_p(wraw.avg.sweepno{1}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
% 
% subplot(2,2,3)
% mapvar=squeeze(mean(wraw.avg.pom{1},1));
% surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
% set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
% set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
% caxset1=caxis;
% caxset1(2)=max(max(mapvar))*0.9;
% caxis(caxset1)
% title (['pom, ' num2str(wraw.avg.sweepno{1}(1)) ' epochs, sig:' num2str(rayleigh_p(wraw.avg.sweepno{1}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
% 
% subplot(2,2,4)
% mapvar=squeeze(mean(wraw.avg.itcm{1},1));
% surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
% set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
% set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
% caxset1=caxis;
% caxset1(2)=max(max(mapvar))*0.9;
% caxis(caxset1)
% title (['itcm, ' num2str(wraw.avg.sweepno{1}(1)) ' epochs, sig:' num2str(rayleigh_p(wraw.avg.sweepno{1}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])

% print ('-djpeg', '-r300', [directory4 'w11_' filenamesout{filecik}(1:end-4) '_' num2str(trigtype(i),'%02.0f') '.jpg']);