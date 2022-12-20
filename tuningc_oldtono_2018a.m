

directory = 'D:\temp\';
% directory = 'D:\Dropbox_Peter\Dropbox\zzz_no_synch\temp\';
filenames = {};

if isempty(filenames)
    filenames=[];
    f=dir( [ directory '*@os.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
end

params_tc.filtere             = [0.5 250];
params_tc.filteru             = [300 6000];
params_tc.high_gamma_filter   = [110 220];
params_tc.filtertype          = 1;
params_tc.newadrate           = 1000;
params_tc.epoch_tframe        = [-25 150];
params_tc.image_tframe        = [-25 150];
params_tc.baseline                 = [-25 0];
params_tc.seleegs.time             = [4 30];

not_tonofile = [];
a_not_tonofile = 0;

bool.artefact_reject    = 1;
bool.eeg                = 0;
bool.convert_sptono     = 0;
params_tc.tonofrq(1)=125;
a = 1;
while params_tc.tonofrq(a)<32000
    a=a+1;
    params_tc.tonofrq(a)=params_tc.tonofrq(a-1)*(2^((1/12)*1));
end

params_tc.calibfrq = params_tc.tonofrq(1:6:end);

params_tc.sptono_limits = [1	7
    8	13
    14	19
    20	25
    26	31
    32	37
    38	43
    44	49
    50	55
    56	61
    62	67
    68	73
    74	79
    80	85
    86	91
    92	97
    ];

% colors1=[0	0	0.5
%     0	0	1
%     0	0.5	1
%     0	0.5	0.5
%     0	0.5	0
%     0	1	0
%     0.5	1	0
%     0.9	0.8	0
%     1	0.5	0
%     0.5	0.5	0
%     0.5	0	0
%     0.5	0	0.5
%     0.5	0	1
%     1	0	1
%     1	0	0];

params_tc.colors2=[0	0	0.5
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
    1	0	0.5
    1   0   1
    1	0.5	1
    0.5	0.5	1
    0.5	0.5 0.5
    0.5	1 0.5
    1	0.5 0.5
    0.1	0.55 0.1
    0.1	0.1	0.55
    0.55	0.1	0.1];

for i1=1:5
params_tc.colors2 = [params_tc.colors2; params_tc.colors2];
end

params_tc.linecolor1=[0.3 0.3 0.3];

params_tc.gridcolor = [0.5 0.5 0.5];

for filecik = 1:length(filenames);
    
    load([directory filenames{filecik}])
    
    [trig_new.trig1s,trig_new.ttype1s, trig_new.triglength1s, findex2s] = module_trig01(trig, params);
    
    trigcik = 1;
    
    trig_new.trig1           = trig_new.trig1s{trigcik};
    trig_new.ttype1          = trig_new.ttype1s{trigcik};
    trig_new.triglength1     = trig_new.triglength1s{trigcik};
    findex2         = findex2s{trigcik};
    
    trig_new.trig1orig       = trig_new.trig1;
    trig_new.trig1           = round(trig_new.trig1./(trig.adrate/params_tc.newadrate));
    
   
    
    if size(trig_new.ttype1,2)>size(trig_new.ttype1,1)
        trig_new.ttype1=trig_new.ttype1';
    end
    
    %%%%%% convert sptono ttype1 to oldtono ttype1
    if max(trig_new.ttype1)>100
        trig_new.ttype2 = trig_new.ttype1;
        trig_new.ttype3 = trig_new.ttype1;
        z = find(trig_new.ttype2>=1000);
        trig_new.ttype2(z)=trig_new.ttype2(z)-1000;
        trig_new.ttype1=trig_new.ttype2;
        if bool.convert_sptono==1
            for i1=1:size(params_tc.sptono_limits,1)
                z = find(trig_new.ttype2>=params_tc.sptono_limits(i1,1) & trig_new.ttype2<=params_tc.sptono_limits(i1,2));
                trig_new.ttype1(z)=i1;
                z = find(trig_new.ttype2>=params_tc.sptono_limits(i1,1)+100 & trig_new.ttype2<=params_tc.sptono_limits(i1,2)+100);
                trig_new.ttype1(z)=i1;
                z = find(trig_new.ttype2>=params_tc.sptono_limits(i1,1)+200 & trig_new.ttype2<=params_tc.sptono_limits(i1,2)+200);
                trig_new.ttype1(z)=i1;
                z = find(trig_new.ttype2>=params_tc.sptono_limits(i1,1)+300 & trig_new.ttype2<=params_tc.sptono_limits(i1,2)+300);
                trig_new.ttype1(z)=i1;
            end
        else
            z = find(trig_new.ttype2>100 & trig_new.ttype2<200);
            trig_new.ttype1(z)=trig_new.ttype1(z)-100;
            z = find(trig_new.ttype2>200 & trig_new.ttype2<300);
            trig_new.ttype1(z)=trig_new.ttype1(z)-200;
            z = find(trig_new.ttype2>300 & trig_new.ttype2<400);
            trig_new.ttype1(z)=trig_new.ttype1(z)-300;
        end
        z = find(trig_new.ttype2==999 | trig_new.ttype2==0 | trig_new.ttype2==1000 | trig_new.ttype2==1999);
        trig_new.ttype1(z)=max(trig_new.ttype1(trig_new.ttype1<100))+1;
        
    else
        trig_new.ttype2 = [];
        trig_new.ttype3 = [];
    end
    
    a=0; trigtype=[];
    for i=0:3000
        z=find(trig_new.ttype1==i);
        if ~isempty(z)
            a=a+1;
            trigtype(a)=i;
        end
    end
    
    if max(trigtype)==15
        trigtype(trigtype==0)=[];
    end
    
    if length(trigtype)>12
        
        params_tc.orig_adrate = craw.adrate;
        adrate  = params_tc.newadrate;
        dt      = 1000/adrate;
        
        [~, ~, cntm, cntc, ~, ~] = module_cnt05(craw, params_tc.newadrate, params_tc.filtere, params_tc.filteru, params_tc.filtertype);
        
        cntc_h      = zeros(size(cntc));
        cntc_hgh    = zeros(size(cntc));
        n = 2;
        Wn = params_tc.high_gamma_filter/(adrate/2);
        [b,a] = butter(n,Wn);
        for chancik=1:size(cntc,1)
            cntc_hgh(chancik,:)=abs(hilbert(filtfilt(b,a,cntc(chancik,:))));
            cntc_h(chancik,:)=abs(hilbert(cntc(chancik,:)));
        end
        
        
        x1 = ceil(params_tc.epoch_tframe(1)/dt);
        x2 = floor(params_tc.epoch_tframe(2)/dt);
        z=find(trig_new.trig1<abs(x1) | trig_new.trig1>length(cntc)-x2);
        trig_new.trig1(z)            = [];
        trig_new.trig1orig(z)        = [];
        trig_new.ttype1(z)           = [];
        trig_new.triglength1(z)      = [];
        
        params_tc.reprate=adrate/median(diff(trig_new.trig1));
        
         eegs.time            = dt*round(params_tc.epoch_tframe(1)/dt):dt:dt*round(params_tc.epoch_tframe(2)/dt);
         
      
             eegs.eegm=zeros(size(cntm,1),length(trig_new.trig1),length(eegs.time));
             x1 = round(params_tc.epoch_tframe(1)/dt);
             x2 = round(params_tc.epoch_tframe(2)/dt);
             for i=1:length(trig_new.trig1)
                 eegs.eegm(:,i,:)=cntm(:,trig_new.trig1(i)+x1:trig_new.trig1(i)+x2);
             end
             
             eegs.eegc=zeros(size(cntc,1),length(trig_new.trig1),length(eegs.time));
             x1 = round(params_tc.epoch_tframe(1)/dt);
             x2 = round(params_tc.epoch_tframe(2)/dt);
             for i=1:length(trig_new.trig1)
                 eegs.eegc(:,i,:)=cntc(:,trig_new.trig1(i)+x1:trig_new.trig1(i)+x2);
             end
             
             eegs.eegc_h=zeros(size(cntc_h,1),length(trig_new.trig1),length(eegs.time));
             x1 = round(params_tc.epoch_tframe(1)/dt);
             x2 = round(params_tc.epoch_tframe(2)/dt);
             for i=1:length(trig_new.trig1)
                 eegs.eegc_h(:,i,:)=cntc_h(:,trig_new.trig1(i)+x1:trig_new.trig1(i)+x2);
             end
             
             eegs.eegc_hgh=zeros(size(cntc_hgh,1),length(trig_new.trig1),length(eegs.time));
             x1 = round(params_tc.epoch_tframe(1)/dt);
             x2 = round(params_tc.epoch_tframe(2)/dt);
             for i=1:length(trig_new.trig1)
                 eegs.eegc_hgh(:,i,:)=cntc_hgh(:,trig_new.trig1(i)+x1:trig_new.trig1(i)+x2);
             end

          
         
         
        
        
        
        if bool.artefact_reject == 1
            %%%%%% artefact reject
            mmax1 = [];
            for i1=1:size(eegs.eegm,2)
                mmax1(i1)=max(mean(abs(squeeze(eegs.eegm(:,i1,:))),1));
            end
            [b,params_tc.deleted_trials1,outliers] = deleteoutliers(mmax1,0.05);
            mmax2 = [];
            for i1=1:size(eegs.eegc,2)
                mmax2(i1)=max(mean(abs(squeeze(eegs.eegc(:,i1,:))),1));
            end
            [b,params_tc.deleted_trials2,outliers] = deleteoutliers(mmax2,0.05);
            mmax3 = [];
            for i1=1:size(eegs.eegc_h,2)
                mmax3(i1)=max(mean(abs(squeeze(eegs.eegc_h(:,i1,:))),1));
            end
            [b,params_tc.deleted_trials3,outliers] = deleteoutliers(mmax3,0.05);
            mmax4 = [];
            for i1=1:size(eegs.eegc_hgh,2)
                mmax4(i1)=max(mean(abs(squeeze(eegs.eegc_hgh(:,i1,:))),1));
            end
            [b,params_tc.deleted_trials4,outliers] = deleteoutliers(mmax4,0.05);
            
            params_tc.deleted_trials = unique([params_tc.deleted_trials1 params_tc.deleted_trials2 params_tc.deleted_trials3 params_tc.deleted_trials4]);
            
            eegs.eegm(:,params_tc.deleted_trials,:)=[];
            eegs.eegc(:,params_tc.deleted_trials,:)=[];
            eegs.eegc_h(:,params_tc.deleted_trials,:)=[];
            eegs.eegc_hgh(:,params_tc.deleted_trials,:)=[];
            trig_new.trig1(params_tc.deleted_trials)            = [];
            trig_new.trig1orig(params_tc.deleted_trials)        = [];
            trig_new.ttype1(params_tc.deleted_trials)           = [];
            trig_new.triglength1(params_tc.deleted_trials)      = [];
        end
        
        % params_tc.baseline
        for i=1:size(eegs.eegm,1)
            for iii=1:size(eegs.eegm,2)
                eegs.eegm(i,iii,:)   = squeeze(eegs.eegm(i,iii,:))-squeeze(mean(eegs.eegm(i,iii,max(find(eegs.time<=params_tc.baseline(1))):max(find(eegs.time<=params_tc.baseline(2)))),3));
                eegs.eegc(i,iii,:)   = squeeze(eegs.eegc(i,iii,:))-squeeze(mean(eegs.eegc(i,iii,max(find(eegs.time<=params_tc.baseline(1))):max(find(eegs.time<=params_tc.baseline(2)))),3));
                eegs.eegc_h(i,iii,:)   = squeeze(eegs.eegc_h(i,iii,:))-squeeze(mean(eegs.eegc_h(i,iii,max(find(eegs.time<=params_tc.baseline(1))):max(find(eegs.time<=params_tc.baseline(2)))),3));
                eegs.eegc_hgh(i,iii,:)   = squeeze(eegs.eegc_hgh(i,iii,:))-squeeze(mean(eegs.eegc_hgh(i,iii,max(find(eegs.time<=params_tc.baseline(1))):max(find(eegs.time<=params_tc.baseline(2)))),3));
            end
        end
        
        
        
        amps.m                  = squeeze(mean(eegs.eegm(:,:,max(find(eegs.time<=params_tc.seleegs.time(1))):max(find(eegs.time<=params_tc.seleegs.time(2)))),3));
        amps.c_h                = squeeze(mean(eegs.eegc_h(:,:,max(find(eegs.time<=params_tc.seleegs.time(1))):max(find(eegs.time<=params_tc.seleegs.time(2)))),3));
        amps.c_hgh              = squeeze(mean(eegs.eegc_hgh(:,:,max(find(eegs.time<=params_tc.seleegs.time(1))):max(find(eegs.time<=params_tc.seleegs.time(2)))),3));
       
        
        trialnumbers = [];
        tunamp = [];
        avgs=[];
        FanoFactor=[];
        
        for trigcik=1:length(trigtype)
            z = find(trig_new.ttype1==trigtype(trigcik));
            trialnumbers(trigcik)=length(z);
            for chancik = 1:size(eegs.eegm,1)
                tunamp.mean.m(chancik,trigcik)    = mean(amps.m(chancik,z),2);
                tunamp.meanp.m(chancik,trigcik)    = signrank(amps.m(chancik,z),zeros(size(amps.m(chancik,z))));
                tunamp.stderr.m(chancik,trigcik)  = std(amps.m(chancik,z),0,2)/sqrt(length(z));
                avgs.m(chancik,trigcik,:)=squeeze(mean(mean(eegs.eegm(chancik,z,:),1),2));
                FanoFactor.m(chancik,trigcik) = (std(amps.m(chancik,z))^2)/mean(amps.m(chancik,z));
                
                tunamp.mean.c_h(chancik,trigcik)    = mean(amps.c_h(chancik,z),2);
                tunamp.meanp.c_h(chancik,trigcik)    = signrank(amps.c_h(chancik,z),zeros(size(amps.c_h(chancik,z))));
                tunamp.stderr.c_h(chancik,trigcik)  = std(amps.c_h(chancik,z),0,2)/sqrt(length(z));
                avgs.c_h(chancik,trigcik,:)=squeeze(mean(mean(eegs.eegc_h(chancik,z,:),1),2));
                FanoFactor.c_h(chancik,trigcik) = (std(amps.c_h(chancik,z))^2)/mean(amps.c_h(chancik,z));
                
                tunamp.mean.c_hgh(chancik,trigcik)    = mean(amps.c_hgh(chancik,z),2);
                tunamp.meanp.c_hgh(chancik,trigcik)    = signrank(amps.c_hgh(chancik,z),zeros(size(amps.c_hgh(chancik,z))));
                tunamp.stderr.c_hgh(chancik,trigcik)  = std(amps.c_hgh(chancik,z),0,2)/sqrt(length(z));
                avgs.c_hgh(chancik,trigcik,:)=squeeze(mean(mean(eegs.eegc_hgh(chancik,z,:),1),2));
                FanoFactor.c_hgh(chancik,trigcik) = (std(amps.c_hgh(chancik,z))^2)/mean(amps.c_hgh(chancik,z));
            end
        end
        avgs.time = eegs.time;
        
        number_of_ptones = max(trigtype)-1;
        
        
        for chancik = 1:size(eegs.eegm,1)
            x = tunamp.mean.m(chancik,1:end-1);
            z=find(x==max(x));
            z=z(1);
            xx =zeros(801,1);
            xx(401-(z-1):401+(number_of_ptones-z))=x;
            tunamp.centeredmean.m(chancik,:)=xx;
            
            x = tunamp.mean.c_h(chancik,1:end-1);
            z=find(x==max(x));
            z=z(1);
            xx =zeros(801,1);
            xx(401-(z-1):401+(number_of_ptones-z))=x;
            tunamp.centeredmean.c_h(chancik,:)=xx;
            
            x = tunamp.mean.c_hgh(chancik,1:end-1);
            z=find(x==max(x));
            z=z(1);
            xx =zeros(801,1);
            xx(401-(z-1):401+(number_of_ptones-z))=x;
            tunamp.centeredmean.c_hgh(chancik,:)=xx;
        end
        tunamp.centeredmean.relativefrq=-400:1:400;
        
        x = tunamp.mean.m(:,1:end-1); %%% to avoid taking noise response into account
        [C,I] = max(x(:));
        [tunamp.max_loc.m(1) tunamp.max_loc.m(2)]=ind2sub(size(x),I);
        tunamp.max.m = C;
        
        x = tunamp.mean.c_h(:,1:end-1); %%% to avoid taking noise response into account
        [C,I] = max(x(:));
        [tunamp.max_loc.c_h(1) tunamp.max_loc.c_h(2)]=ind2sub(size(x),I);
        tunamp.max.c_h = C;
        
        x = tunamp.mean.c_hgh(:,1:end-1); %%% to avoid taking noise response into account
        [C,I] = max(x(:));
        [tunamp.max_loc.c_hgh(1) tunamp.max_loc.c_hgh(2)]=ind2sub(size(x),I);
        tunamp.max.c_hgh = C;
        
        tunamp_mean = tunamp.mean.m;
        FF = FanoFactor.m;
        for i1=1:size(tunamp_mean,1)
            [M,I]=max(tunamp_mean(i1,1:end-1));
            lam_FF(i1)=FF(i1,I);
        end
        FanoFactor.lam.m = lam_FF;
        
        I=[]; lam_FF =[];
        tunamp_mean = tunamp.mean.m;
        FF = FanoFactor.m;
        for i1=1:size(tunamp_mean,1)
            [M,I(i1)]=max(tunamp_mean(i1,1:end-1));
            lam_FF(i1)=FF(i1,I(i1));
        end
        FanoFactor.lam.m = lam_FF;
        FanoFactor.lam.mloc = I;
        
        I=[]; lam_FF =[];
        tunamp_mean = tunamp.mean.c_h;
        FF = FanoFactor.c_h;
        for i1=1:size(tunamp_mean,1)
            [M,I(i1)]=max(tunamp_mean(i1,1:end-1));
            lam_FF(i1)=FF(i1,I(i1));
        end
        FanoFactor.lam.c_h = lam_FF;
        FanoFactor.lam.c_hloc = I;
        
        I=[]; lam_FF =[];
        tunamp_mean = tunamp.mean.c_hgh;
        FF = FanoFactor.c_hgh;
        for i1=1:size(tunamp_mean,1)
            [M,I(i1)]=max(tunamp_mean(i1,1:end-1));
            lam_FF(i1)=FF(i1,I(i1));
        end
        FanoFactor.lam.c_hgh = lam_FF;
        FanoFactor.lam.c_hghloc = I;
        
        curfig = figure;
        set(curfig,'position',[100   60   1272   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
        figureax1=axes('Position',[0 0 1 1],'Visible','off');
        set(curfig,'CurrentAxes',figureax1)
        text(0.5,0.985,[filenames{filecik}],'FontSize',10,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
        
        fsize = 6;
        
        axes('position',[0.04 0.89 0.14 0.07],'box','off')
        tunamp_mean = tunamp.mean.m(tunamp.max_loc.m(1),:);
        tunamp_stderr = tunamp.stderr.m(tunamp.max_loc.m(1),:);
        for i=1:size(tunamp_mean,2)
            errorbar(i,tunamp_mean(i),tunamp_stderr(i),'color',params_tc.colors2(i,:),'linewidth',2)
            hold on
        end
        plot(1:size(tunamp_mean,2),tunamp_mean,'color',params_tc.linecolor1)
        set(gca,'xlim',[0.5 size(tunamp_mean,2)+0.5])
        set(gca,'xtick',1:1:size(tunamp_mean,2),'xticklabel',trigtype,'box','off')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['MUA ch ' num2str(tunamp.max_loc.m(1)) ', frq ' num2str(tunamp.max_loc.m(2)) ', FF ' num2str(FanoFactor.m(tunamp.max_loc.m(1),tunamp.max_loc.m(2)))])
        
        axes('position',[0.19 0.89 0.14 0.07],'box','off')
        x=squeeze(avgs.m(tunamp.max_loc.m(1),:,max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))));
        ylim2=max(max(x))*1.05;
        if min(min(x))>0
            ylim1=min(min(x))*0.95;
        else
            ylim1=min(min(x))*1.05;
        end
        set(gca,'NextPlot','replacechildren','colororder',params_tc.colors2)
        plot(eegs.time(max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))),x)
        set(gca,'xlim',[params_tc.image_tframe(1) params_tc.image_tframe(2)],'ylim',[ylim1 ylim2])
        line([params_tc.seleegs.time(1) params_tc.seleegs.time(1)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        line([params_tc.seleegs.time(2) params_tc.seleegs.time(2)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['eegs.time interval: ' num2str(params_tc.seleegs.time(1)) ' - ' num2str(params_tc.seleegs.time(2))])
        
        axes('position',[0.36 0.89 0.14 0.07],'box','off')
        tunamp_mean = tunamp.mean.c_h(tunamp.max_loc.c_h(1),:);
        tunamp_stderr = tunamp.stderr.c_h(tunamp.max_loc.c_h(1),:);
        for i=1:size(tunamp_mean,2)
            errorbar(i,tunamp_mean(i),tunamp_stderr(i),'color',params_tc.colors2(i,:),'linewidth',2)
            hold on
        end
        plot(1:size(tunamp_mean,2),tunamp_mean,'color',params_tc.linecolor1)
        set(gca,'xlim',[0.5 size(tunamp_mean,2)+0.5])
        set(gca,'xtick',1:1:size(tunamp_mean,2),'xticklabel',trigtype,'box','off')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['CSD ch ' num2str(tunamp.max_loc.c_h(1)) ', frq ' num2str(tunamp.max_loc.c_h(2)) ', FF ' num2str(FanoFactor.c_h(tunamp.max_loc.c_h(1),tunamp.max_loc.c_h(2)))])
        
        axes('position',[0.51 0.89 0.14 0.07],'box','off')
        x=squeeze(avgs.c_h(tunamp.max_loc.c_h(1),:,max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))));
        ylim2=max(max(x))*1.05;
        if min(min(x))>0
            ylim1=min(min(x))*0.95;
        else
            ylim1=min(min(x))*1.05;
        end
        set(gca,'NextPlot','replacechildren','colororder',params_tc.colors2)
        plot(eegs.time(max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))),x)
        set(gca,'xlim',[params_tc.image_tframe(1) params_tc.image_tframe(2)],'ylim',[ylim1 ylim2])
        line([params_tc.seleegs.time(1) params_tc.seleegs.time(1)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        line([params_tc.seleegs.time(2) params_tc.seleegs.time(2)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['eegs.time interval: ' num2str(params_tc.seleegs.time(1)) ' - ' num2str(params_tc.seleegs.time(2))])
        
        axes('position',[0.68 0.89 0.14 0.07],'box','off')
        tunamp_mean = tunamp.mean.c_hgh(tunamp.max_loc.c_hgh(1),:);
        tunamp_stderr = tunamp.stderr.c_hgh(tunamp.max_loc.c_hgh(1),:);
        for i=1:size(tunamp_mean,2)
            errorbar(i,tunamp_mean(i),tunamp_stderr(i),'color',params_tc.colors2(i,:),'linewidth',2)
            hold on
        end
        plot(1:size(tunamp_mean,2),tunamp_mean,'color',params_tc.linecolor1)
        set(gca,'xlim',[0.5 size(tunamp_mean,2)+0.5])
        set(gca,'xtick',1:1:size(tunamp_mean,2),'xticklabel',trigtype,'box','off')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['HighGamma ch ' num2str(tunamp.max_loc.c_hgh(1)) ', frq ' num2str(tunamp.max_loc.c_hgh(2)) ', FF ' num2str(FanoFactor.c_hgh(tunamp.max_loc.c_hgh(1),tunamp.max_loc.c_hgh(2)))])
        
        axes('position',[0.83 0.89 0.14 0.07],'box','off')
        x=squeeze(avgs.c_hgh(tunamp.max_loc.c_hgh(1),:,max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))));
        ylim2=max(max(x))*1.05;
        if min(min(x))>0
            ylim1=min(min(x))*0.95;
        else
            ylim1=min(min(x))*1.05;
        end
        set(gca,'NextPlot','replacechildren','colororder',params_tc.colors2)
        plot(eegs.time(max(find(eegs.time<=params_tc.image_tframe(1))):max(find(eegs.time<=params_tc.image_tframe(2)))),x)
        set(gca,'xlim',[params_tc.image_tframe(1) params_tc.image_tframe(2)],'ylim',[ylim1 ylim2])
        line([params_tc.seleegs.time(1) params_tc.seleegs.time(1)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        line([params_tc.seleegs.time(2) params_tc.seleegs.time(2)],get(gca,'ylim'),'color',params_tc.linecolor1,'linestyle',':')
        set(gca,'xcolor',[.1 .1 .1],'ycolor',[.1 .1 .1],'color',[1 1 1],'fontsize',fsize)
        title (['eegs.time interval: ' num2str(params_tc.seleegs.time(1)) ' - ' num2str(params_tc.seleegs.time(2))])
        
        axpos=[0.04 0.04 0.14 0.82];
        z = find (trig_new.ttype1 ==  tunamp.max_loc.m(2));
        figureax2a=axes('Position',axpos);
        [cax1] = csd_maker_no_subplot07(squeeze(mean(eegs.eegc(:,z,:),2)),eegs.time,1,params_tc.image_tframe,[0 0],[],axpos,figureax2a);
        % line(eraw.eegs.timeframe,[selchan selchan],'LineStyle',':','LineWidth',0.5,'Color',[0 0 .8])
        title([ 'CSD ,   cax: ' num2str(cax1(1)) ' - ' num2str(cax1(2))],'FontSize',fsize)
        colormap(figureax2a,flipud(jet))
        
        axpos=[0.22 0.04 0.22 0.82];
        figureax6a=axes('Position',axpos);
        stepm = 1;
        tunn = [];
        yticks = [];
        tunmax = [];
        tun = tunamp.mean.m(:,:);
        for i1=1:size(tun,1)
            tunn(i1,:)=tun(i1,:)-min(tun(i1,:));
            tunn(i1,:)=tunn(i1,:).*1/max(tunn(i1,:));
            x=find(tunn(i1,:)==max(tunn(i1,:)));
            tunmax(i1)=x(1);
        end
        
        
        for i1=1:size(tun,1)
            plot(1:size(tunn,2),tunn(i1,:)-stepm*(i1-1),'color',params_tc.colors2(i1,:),'linewidth',1)
            set(gca,'xlim',[0.5 size(tunn,2)+0.5])
            yticks(i1)=-stepm*(i1-1);
            hold on
            line(get(gca,'xlim'),[yticks(i1) yticks(i1)],'color',params_tc.gridcolor,'linestyle',':')
        end
        set(gca,'ylim',[-20 1])
        ylabels     = num2str(tunmax','%11.3g');
        clear ylabels2
        for i1=1:size(ylabels,1)
            ylabels2(i1,:)=[num2str(i1,'%02.0f') '  ' ylabels(i1,:)];
        end
        ylabels2=flipud(ylabels2);
        set(gca,'ytick',fliplr(yticks),'yticklabel',ylabels2,'fontsize',fsize,'xtick',1:round(size(tunn,2)/15):size(tunn,2))
        title(['MUA tuning'],'FontSize',fsize)
        
        axpos=[0.48 0.04 0.22 0.82];
        figureax6a=axes('Position',axpos);
        stepm = 1;
        tunn = [];
        yticks = [];
        tunmax = [];
        tun = tunamp.mean.c_h(:,:);
        for i1=1:size(tun,1)
            tunn(i1,:)=tun(i1,:)-min(tun(i1,:));
            tunn(i1,:)=tunn(i1,:).*1/max(tunn(i1,:));
            x=find(tunn(i1,:)==max(tunn(i1,:)));
            tunmax(i1)=x(1);
        end
        
        
        for i1=1:size(tun,1)
            plot(1:size(tunn,2),tunn(i1,:)-stepm*(i1-1),'color',params_tc.colors2(i1,:),'linewidth',1)
            set(gca,'xlim',[0.5 size(tunn,2)+0.5])
            yticks(i1)=-stepm*(i1-1);
            hold on
            line(get(gca,'xlim'),[yticks(i1) yticks(i1)],'color',params_tc.gridcolor,'linestyle',':')
        end
        set(gca,'ylim',[-20 1])
        ylabels     = num2str(tunmax','%11.3g');
        clear ylabels2
        for i1=1:size(ylabels,1)
            ylabels2(i1,:)=[num2str(i1,'%02.0f') '  ' ylabels(i1,:)];
        end
        ylabels2=flipud(ylabels2);
        set(gca,'ytick',fliplr(yticks),'yticklabel',ylabels2,'fontsize',fsize,'xtick',1:round(size(tunn,2)/15):size(tunn,2))
        title(['CSD tuning'],'FontSize',fsize)
        
        axpos=[0.74 0.04 0.22 0.82];
        figureax6a=axes('Position',axpos);
        stepm = 1;
        tunn = [];
        yticks = [];
        tunmax = [];
        tun = tunamp.mean.c_hgh(:,:);
        for i1=1:size(tun,1)
            tunn(i1,:)=tun(i1,:)-min(tun(i1,:));
            tunn(i1,:)=tunn(i1,:).*1/max(tunn(i1,:));
            x=find(tunn(i1,:)==max(tunn(i1,:)));
            tunmax(i1)=x(1);
        end
        
        for i1=1:size(tun,1)
            plot(1:size(tunn,2),tunn(i1,:)-stepm*(i1-1),'color',params_tc.colors2(i1,:),'linewidth',1)
            set(gca,'xlim',[0.5 size(tunn,2)+0.5])
            yticks(i1)=-stepm*(i1-1);
            hold on
            line(get(gca,'xlim'),[yticks(i1) yticks(i1)],'color',params_tc.gridcolor,'linestyle',':')
        end
        set(gca,'ylim',[-20 1])
        ylabels     = num2str(tunmax','%11.3g');
        clear ylabels2
        for i1=1:size(ylabels,1)
            ylabels2(i1,:)=[num2str(i1,'%02.0f') '  ' ylabels(i1,:)];
        end
        ylabels2=flipud(ylabels2);
        set(gca,'ytick',fliplr(yticks),'yticklabel',ylabels2,'fontsize',fsize,'xtick',1:round(size(tunn,2)/15):size(tunn,2))
        title(['HighGamma tuning'],'FontSize',fsize)
        
        if bool.artefact_reject == 1 & bool.convert_sptono == 1
            
            xx = 'a';
        elseif bool.convert_sptono == 0
             xx = 'b';
        else
            xx = [];
        end
        
        print ('-djpeg', '-r300', [directory 'tun01' xx '_' filenames{filecik}(1:end-7) '.jpg']);
        
        close all
        
        clear mmax1 mmax2 mmax3 mmax4 a b axpos C cax1 chancik clf craw cai cntc cntc_h cntc_hgh cntm f ff i i1 iii n outliers stepm tunamp_mean tunamp_stderr trig I dt adrate fsize tun tunn x x1 x2 ylim1 ylim2 z
        
        if bool.eeg < 1
             eegs.eegm      = [];
             eegs.eegc      = [];
             eegs.eegc_h    = [];
             eegs.eegc_hgh  = [];
        end
        
        save ([directory 'tun01' xx '_' filenames{filecik}(1:end-7) '.mat'],'amps','avgs','bool','eegs','FanoFactor','params_tc','trig_new','tunamp','params' ,'-v7.3');
        
    else
        a_not_tonofile = a_not_tonofile+1;
        not_tonofile(a_not_tonofile)=filecik;
    end
end


