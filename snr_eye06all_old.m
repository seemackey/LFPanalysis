
function [eyedata] = snr_eye06all_old(directory,det_met)

%%%% examples:

% directory = 'E:\dyneyen2\bbn\BBN_eo\';
% det_met   = [2 0];
% [eyedata] = snr_eye06all_old(directory, det_met);

filenames = {};

tic

if isempty(filenames)
    filenames=[];
    f=dir( [ directory '*@os*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
end

%%%%%% parameters for all methods
bool.detection_method   = 6;    % 6 is just the newest version

bool.merge_window   = 50;
bool.fixdur         = 100;                        %%% minimum fixation duration (ms)
bool.merge          = 1;

%%%%%%%%%% switch on and off the different methods
if isempty(det_met)
    bool.det_met(1) = 2; % 33
    bool.det_met(2) = 2;  % 2
else
    bool.det_met=det_met;
end



bool.min_segment_length = 1000;         %%% minimum length of a segment without eyes closed or blink periods (ms)
bool.newadrate          = 100;         %%% sampling rate for saccade data (Hz)
bool.eyesclosed_edge    = 100;            %%% edges around blinks and closed eyes (ms)
bool.eyesclosed_limit   = 1800;         %%% limit for bad signal (amp)
bool.isi_limit          = 1250;         %%% limit for exluding saccades from pooled isi (ms)

%%%%%% parameters for method 2 (Engbert-Kliegl algorithm)
bool.parameter_method2                  = 3.2;     %%% 3.2 seems to work for most files
bool.sacc_duration_method2              = [2 8];    %%% 8 seems to be ok
bool.method2arej                        = 1;

%%%%% parameters for method3  (ClusterFix)
bool.duration_threshold_for_saccades    = 10;     %%% minimum saccade duration (ms) default was 8, now trying 10
bool.duration_threshold_for_fixations   = 50;      %%% minimum fixation duration (ms)
bool.fltord                             = 60;
bool.lowpasfrq                          = 100;

%%%%%%%%%%%%%%%%% artefact reject (cleanup) switches

bool.cleanup4           = 1; %%%%%% based on sacc amp, and pre and post sacc amp
bool.saccamp_ratios     = 0.8; %%%%%% 0.5 for new data, 0.8 for old
%%%%%%%%%%%%%%%% figure and data saving options

bool.file_edge  = 250;
bool.image      = 2;     % 1 if show, 2 if close, 0 if do not even show
bool.print      = 1;     % if 1, it will save the image as jpeg
bool.inject     = 2;     % add led eyes open and saccade start, end triggers in anatrig positions 3, 4 and 5 respectively
% 1 just "injects the file, 2 saves a new instance
% renamed based on method
bool.saveraw    = 1;     % 1 if all eyedata at 1000 Hz is needed

waitb           = waitbar(0);

for filecik = 1:length(filenames)
    % for filecik = 1:1
    
    a2=['file' num2str(filecik) ' - ' filenames{filecik}];
    waitbar(filecik/length(filenames) , waitb,a2);
    
    load([directory filenames{filecik}],'trig','cai','params')
    
    
    
    cnts.eyepos         = cai.cnt(2:4,:); %%% get the eye data
    
    nyq = bool.newadrate/2.01;
    n = 2;
    Wn = nyq/(cai.adrate/2);
    [b,a] = butter(n,Wn,'low');
    for i=1:size(cnts.eyepos,1)
        cnts.eyepos(i,:)=filtfilt(b,a,cnts.eyepos(i,:));
        % cnts.eyepos(i,:)=smoothdata(cnts.eyepos(i,:),'gaussian',20);
    end
    cnts.eyepos      = resample(cnts.eyepos',bool.newadrate,cai.adrate)';
    adrate_orig     = cai.adrate;
    cnts.adrate      = bool.newadrate;
    eyedata.adrate   = bool.newadrate;
    
    %%%%% original closed criterion, x or y lost is considered closed
    eyesclosed    = find(abs(cnts.eyepos(1,:))>bool.eyesclosed_limit | abs(cnts.eyepos(2,:))>bool.eyesclosed_limit);
    
    %%%%% stricter eyes closed criterion, both x and y have to be lost for eyes closed
    % eyesclosed      = find(abs(cnts.eyepos(1,:))>bool.eyesclosed_limit & abs(cnts.eyepos(2,:))>bool.eyesclosed_limit);
    
    edge1           = round(cnts.adrate*bool.eyesclosed_edge/1000);
    [seg_eyesclosed, segvec_eyesclosed]   = windows02(eyesclosed',1,edge1,length(cnts.eyepos(1,:)));
    eyepos = cnts.eyepos;
    eyepos(:,segvec_eyesclosed)=NaN;
    
    %%%% NaN the edges
    eyepos(:,1:bool.file_edge)=NaN;
    eyepos(:,end-bool.file_edge:end)=NaN;
    
    eyesopen      = find(~isnan(eyepos(1,:)));
    edge2 = 0;
    min_segment_length = round(cnts.adrate*bool.min_segment_length/1000);
    [seg_eyesopen, segvec_eyesopen]   = windows02(eyesopen',min_segment_length,edge2,length(cnts.eyepos(1,:)));
    
    eyeposdiff = diff(eyepos')';
    
    eyepos(end+1,:)=1:size(cnts.eyepos,2);
    
    angles = [];
    lengths = [];
    for tcik = 1:size(eyeposdiff,2)
        z=complex(eyeposdiff(1,tcik),eyeposdiff(2,tcik));
        angles(tcik) = angle(z);
        lengths(tcik) = abs(z);
    end
    
    angles = [NaN, angles];
    lengths = [NaN, lengths];
    
    eyepos_x=cnts.eyepos(1,:);
    eyepos_y=cnts.eyepos(2,:);
    
    ec = find(isnan(eyepos(1,:)) & isnan(eyepos(2,:)));
    eo = find(~isnan(eyepos(1,:)) & ~isnan(eyepos(2,:)));
    
    eyedata.ec_eo_ratio(filecik) = length(ec)/length(eo);
    
    %%%%%% led related stuff
    eyedata.ledtrig_eo                  = [];
    eyedata.ledtrig_ec                  = [];
    
    trig1_ec = [];
    trig1_eo = [];
    a1=0;
    a2=0;
    %%%%%%%%%%%% this is the part you have to play with replicate this for
    %%%%%%%%%%%% trig.anatrig{1}
    if ~isempty(trig.anatrig{2})
        trig1 = round(trig.anatrig{2}/(trig.adrate/bool.newadrate));
        for i1=1:length(trig1)
            if find(trig1(i1)==eo)
                a1=a1+1;
                trig1_eo(a1)=i1;
            else
                a2=a2+1;
                trig1_ec(a2)=i1;
            end
        end
        eyedata.ledtrig_eo                = trig1_eo;
        eyedata.ledtrig_ec                 = trig1_ec;
        eyedata.ledtrig_ec_eo_ratio = length(eyedata.ledtrig_ec)/length(eyedata.ledtrig_eo);
    end
    
    trig1_ec = [];
    trig1_eo = [];
    a1=0;
    a2=0;
    if ~isempty(trig.anatrig{1})
        trig1 = round(trig.anatrig{1}/(trig.adrate/bool.newadrate));
        for i1=1:length(trig1)
            if find(trig1(i1)==eo)
                a1=a1+1;
                trig1_eo(a1)=i1;
            else
                a2=a2+1;
                trig1_ec(a2)=i1;
            end
        end
        eyedata.audtrig_eo                = trig1_eo;
        eyedata.audtrig_ec                 = trig1_ec;
        eyedata.audtrig_ec_eo_ratio = length(eyedata.ledtrig_ec)/length(eyedata.ledtrig_eo);
        
    end
    
 
    
    methodcik = 0;
    
    if bool.det_met(1) == 2
        
        samprate=1/bool.newadrate;
        
        ss1 = [];
        ss2 = [];
        %%%%%%% method 33
        for segcik = 1:size(seg_eyesopen,1)
            
            sac = [];
            
            seg         = cnts.eyepos(:,seg_eyesopen(segcik,1):seg_eyesopen(segcik,2));
            seg_time    = eyepos(4,seg_eyesopen(segcik,1):seg_eyesopen(segcik,2));
            
            
            
            [fixationstats]     = ClusterFix01({seg},samprate, bool.duration_threshold_for_saccades, bool.duration_threshold_for_fixations,bool.fltord,bool.lowpasfrq);
            
            if ~isempty(fixationstats{1,1}.saccadetimes)
                sac(:,1)=fixationstats{1,1}.saccadetimes(1,:)';
                sac(:,2)=fixationstats{1,1}.saccadetimes(2,:)';
                sac(:,3:end+6)=fixationstats{1,1}.SaaccadeClusterValues(:,:);
            else
                sac=[];
            end
            
            if ~isempty(sac)
                if sac(1,1)==1
                    sac(1,:)=[];
                end
                sac(:,1)=seg_time(sac(:,1));
                sac(:,2)=seg_time(sac(:,2));
            end
            
            ss1     = [ ss1; sac];
        end
        
        methodcik = methodcik+1;
        ss1s{methodcik}=ss1;
    end
    
    
    
    if bool.det_met(2) == 2
        ss1 = [];
        %%%%%%% method 2
        for segcik = 1:size(seg_eyesopen,1)
            
            seg         = cnts.eyepos(1:2,seg_eyesopen(segcik,1):seg_eyesopen(segcik,2));
            seg_time    = eyepos(4,seg_eyesopen(segcik,1):seg_eyesopen(segcik,2));
            
            % x=smoothdata(seg');
            x=seg';
            % x(1,:)=[];
            [sac, radius] = microsacc(x,bool.parameter_method2,bool.sacc_duration_method2(1),bool.newadrate);
            
            if ~isempty(sac)
                if sac(1,1)==1
                    sac(1,:)=[];
                end
                sac(:,1)=seg_time(sac(:,1));
                sac(:,2)=seg_time(sac(:,2));
                
                %%%% this method does not work without this correction
                if bool.method2arej == 1
                    xdiff=[NaN diff(sac(:,1))'];
                    dels1=find(xdiff<bool.sacc_duration_method2(2));
                    sac(dels1,:)=[];
                end
            end
            
            ss1     = [ ss1; sac];
        end
        
        % ss1(:,1:2)= ss1(:,1:2).*1000/eyedata.adrate;
        
        methodcik = methodcik+1;
        ss1s{methodcik}=ss1;
        
    end
    
    if ~isempty(ss1s{1})
    if length(ss1s)>1 & ~isempty(ss1s{2})
    x1 = ss1s{1}(:,1:2);
    x2 = ss1s{2}(:,1:2);
     sacc1=x1(:,1);
    sacc2=x2(:,1);
    else
        x1 = ss1s{1}(:,1:2);
        x2=[];
         sacc1=x1(:,1);
    sacc2=[];
    end
    else
        sacc1=[];
        sacc2=[];
        x1=[];
        x2=[];
    end
    
    
    

    sacctimediff = [];
    zs = [];
    for i1=1:length(sacc2)
        sacctimediff(i1)=min(abs(sacc2(i1)-sacc1));
        zsacc = find(abs(sacc2(i1)-sacc1)==sacctimediff(i1));
        zs(i1)=zsacc(1);
    end
    zsacc = find(sacctimediff<bool.merge_window);
    eyedata.coincidence_ratio=length(zsacc)/length(sacc1);
    
    try
    x2(zsacc,:)=[];
    catch
        x2=[];
    end
    
    saccades0=[x1; x2];
    
    try
    [B,I] = sort(saccades0(:,1));
    saccades0 = saccades0(I,:);
    catch
        saccades0 = saccades0(:,:);
    end
    
    saccades0_origin = [];
    saccades0_origin(:,1) = [ones(length(x1),1); ones(length(x2),1)+1];
    saccades0_origin(:,2) = [1:length(x1), 1:length(x2)];
    
    try
    saccades0_origin = saccades0_origin(I,:);
    catch
        saccades0_origin = saccades0_origin(:,:);
    end
    
    eyedata.saccades0 = saccades0;
    eyedata.saccades0_origin = saccades0_origin;
    if ~isempty(ss1s{1})
    if det_met(1)==0
        eyedata.saccnums = [0 length(ss1s{1}(:,1:2)) length(zsacc) length(eyedata.saccades0)];
    elseif det_met(2)==0
        eyedata.saccnums = [length(ss1s{1}(:,1:2)) 0 length(zsacc) length(eyedata.saccades0)];
    else
        eyedata.saccnums = [length(ss1s{1}(:,1:2)) length(ss1s{2}(:,1:2)) length(zsacc) length(eyedata.saccades0)];
    end
    else
        eyedata.saccnums = [0 0 0 0];
    end
    
    eyedata.saccamps        = [];
    eyedata.saccamps1a        = [];
    eyedata.saccamps1b        = [];
    eyedata.saccamps2a        = [];
    eyedata.saccamps2b        = [];
    eyedata.postsaccamps1    = [];
    eyedata.presaccamps1    = [];
    for i1=1:size(saccades0,1)
        eyedata.saccamps(i1) = nanmean(lengths(saccades0(i1,1):saccades0(i1,2)));
        eyedata.saccamps1a(i1) = nanmean(diff(cnts.eyepos(1,saccades0(i1,1):saccades0(i1,2))));
        eyedata.saccamps1b(i1) = nanmean(diff(cnts.eyepos(2,saccades0(i1,1):saccades0(i1,2))));
        eyedata.saccamps2a(i1) = nanmean(abs(diff(cnts.eyepos(1,saccades0(i1,1):saccades0(i1,2)))));
        eyedata.saccamps2b(i1) = nanmean(abs(diff(cnts.eyepos(2,saccades0(i1,1):saccades0(i1,2)))));
        try
            seg00 = cnts.adrate*bool.fixdur/1000;
            eyedata.postsaccamps1(i1) = nanmean(lengths(saccades0(i1,2)+1:saccades0(i1,2)+seg00));
            eyedata.presaccamps1(i1) = nanmean(lengths(saccades0(i1,1)-seg00-1:saccades0(i1,1)-1));
            
        catch
            eyedata.postsaccamps1(i1) = 0;
            eyedata.presaccamps1(i1) = 0;
        end
    end
    
    
    dels0        = [];
    dels1        = [];
    dels2        = [];
    dels3        = [];
    dels4        = [];
    
    %%%%%% cleanup4, based on saccade velocity
    if ~isempty(eyedata.saccamps)
        if bool.cleanup4 == 1
            dels4 = find(eyedata.postsaccamps1>eyedata.saccamps.*bool.saccamp_ratios | eyedata.presaccamps1>eyedata.saccamps.*bool.saccamp_ratios);
        end
    end
    
    %%%%%%%%% deleted and selected "saccades"
    try
        dels        = sort(unique([dels0; dels1; dels2; dels3; dels4]));
    catch
        dels        = sort(unique([dels0; dels1'; dels2; dels3; dels4]));
    end
    

    if ~isempty(saccades0)

        xx = [NaN; diff(saccades0(:,1))]./bool.newadrate*1000;

    else
        xx = [];
    end
    eyedata.intersacc_interval = xx;
    
    if ~isempty(saccades0)
        eyedata.saccdur = (saccades0(:,2)-saccades0(:,1))./eyedata.adrate*1000;
    else
        eyedata.saccdur = [];
    end
    
    
    ec = find(isnan(eyepos(1,:)));
    eo = find(~isnan(eyepos(1,:)));
    
    eyedata.ec_eo_ratio = length(ec)/length(eo);
    
    eyedata.segs.seg_eyesopen       = seg_eyesopen;
    eyedata.segs.segvec_eyesopen    = segvec_eyesopen;
    
    eyedata.segs.seg_eyesclosed     = seg_eyesclosed;
    eyedata.segs.segvec_eyesclosed  = segvec_eyesclosed;
    
    aa = sprintf('%d', bool.det_met);
    
    fname = [filenames{filecik}(1:end-4) '_' 'eye0' num2str(bool.detection_method) '_' aa];

   
    %%%%% inject file
    if ~isempty(saccades0)
        if bool.inject == 1
            m = matfile([directory filenames{filecik}],'Writable',true);
            
            trig = m.trig;
            if isempty(trig.trigchan)
                craw = m.craw;
                trig.adrate = craw.adrate;
            end
            if ~isempty(trig.anatrig{2})
                trig.anatrig{3}=trig.anatrig{2}(eyedata.ledtrig_eo);
            end
            
            if ~isempty(trig.anatrig{1})
               
                trig.anatrig{3}=trig.anatrig{1}(eyedata.audtrig_eo);
            end
            
            trig.anatrig{4}=eyedata.saccades0(:,1)*(trig.adrate/bool.newadrate);
            trig.anatrig{5}=eyedata.saccades0(:,2)*(trig.adrate/bool.newadrate);
            trig.trigchan = [1:5];
            params = m.params;
            
            params.saccade                      = [];
            
            params.saccade.eyedata              = eyedata;
            params.saccade.bool                 = bool;
            
            %%%% artefact rejection if applicable
            params.saccade.eyedata.arej         = ones(1,size(saccades0,1));
            params.saccade.eyedata.arej(dels)   = 0;
            
            params.saccade.eyedata.ss1s          = ss1s;
            
            
            if bool.saveraw == 1
                params.saccade.eyedata.lengths = lengths;
                params.saccade.eyedata.angles  = angles;
                params.saccade.eyedata.eyepos  = eyepos;
                params.saccade.eyedata.cnt      = cnts;
            end
            
            m.params = params;
            m.trig = trig;
            clear m
        elseif bool.inject == 2
            xx = whos('-file',[directory filenames{filecik}]);
            load([directory filenames{filecik}])
            
            if isempty(trig.trigchan)
                trig.adrate = craw.adrate;
            end
            if ~isempty(trig.anatrig{2})
                trig.anatrig{3}=trig.anatrig{2}(eyedata.ledtrig_eo);
            end
            
            if ~isempty(trig.anatrig{1})
               
                trig.anatrig{3}=trig.anatrig{1}(eyedata.audtrig_eo);
            end
            
            trig.anatrig{4}=eyedata.saccades0(:,1)*(trig.adrate/bool.newadrate);
            trig.anatrig{5}=eyedata.saccades0(:,2)*(trig.adrate/bool.newadrate);
            trig.trigchan = [1:5];
            
            params.saccade                      = [];
            
            params.saccade.eyedata              = eyedata;
            params.saccade.bool                 = bool;
            
            %%%% artefact rejection if applicable
            params.saccade.eyedata.arej         = ones(1,size(saccades0,1));
            params.saccade.eyedata.arej(dels)   = 0;
            
            params.saccade.eyedata.ss1s          = ss1s;
            
            if bool.saveraw == 1
                params.saccade.eyedata.lengths = lengths;
                params.saccade.eyedata.angles  = angles;
                params.saccade.eyedata.eyepos  = eyepos;
                params.saccade.eyedata.cnt      = cnts;
            end
            
            save([directory fname '.mat'],xx.name,'-v7.3')
        else
            
            trig.anatrig{4}=eyedata.saccades0{filecik}(:,1)*(trig.adrate/bool.newadrate);
            trig.anatrig{5}=eyedata.saccades0{filecik}(:,2)*(trig.adrate/bool.newadrate);
            trig.trigchan = [1:5];
        end
        
        
        %%%%%%%%%% plotting
        
        
        
        if bool.image > 0
            
            
            
            curfig = figure;
            set(curfig,'position',[100   50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
            
            if bool.image == 11
                plot(eyedata.eyepos{filecik}(1:2,:)')
                hold on
                
                zp=find(params.saccade.eyedata.arej == 1);
                
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),1) eyedata.saccades0(zp(i1),1)],get(gca,'ylim'),'color',[1 0 0],'linestyle',':')
                end
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),2) eyedata.saccades0(zp(i1),2)],get(gca,'ylim'),'color',[0 1 0],'linestyle',':')
                end
                
                %         if ~isempty(trig.anatrig{2})
                %             x = trig.anatrig{2}(eyedata.ledtrig_eo{filecik})./(trig.adrate/bool.newadrate);
                %             for i1=1:length(x)
                %                 line([x(i1) x(i1)],get(gca,'ylim'),'color',[0 0 .5],'linestyle',':')
                %             end
                %         end
            else
                
                subplot(2,3,1)
                plot(params.saccade.eyedata.lengths)
                hold on
                
                zp=find(params.saccade.eyedata.arej == 1);
                
                try
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),1) eyedata.saccades0(zp(i1),1)],get(gca,'ylim'),'color',[1 0 0],'linestyle',':')
                end
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),2) eyedata.saccades0(zp(i1),2)],get(gca,'ylim'),'color',[0 1 0],'linestyle',':')
                end
                if ~isempty(trig.anatrig{2})
                    x = trig.anatrig{2}(eyedata.ledtrig_eo)./(trig.adrate/bool.newadrate);
                    for i1=1:length(x)
                        line([x(i1) x(i1)],get(gca,'ylim'),'color',[0 0 .5],'linestyle',':')
                    end
                end
                catch
                end
                title(num2str(eyedata.saccnums))
                subplot(2,3,2)
                plot(params.saccade.eyedata.eyepos(1:2,:)')
                hold on
                zp=find(params.saccade.eyedata.arej == 1);
                
                try
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),1) eyedata.saccades0(zp(i1),1)],get(gca,'ylim'),'color',[1 0 0],'linestyle',':')
                end
                for i1=1:length(zp)
                    line([eyedata.saccades0(zp(i1),2) eyedata.saccades0(zp(i1),2)],get(gca,'ylim'),'color',[0 1 0],'linestyle',':')
                end
                catch
                end
        
                subplot(2,3,4)
                hist(eyedata.saccdur,0:1:150)
                set(gca,'xlim',[0 150])
                title(['median saccdur: ' num2str(median(eyedata.saccdur))])
                
                subplot(2,3,5)
                hist(eyedata.intersacc_interval,0:50:2000)
                set(gca,'xlim',[0 1000])
                z = find(eyedata.intersacc_interval<=bool.isi_limit);
                title(['median isi ' num2str(median(eyedata.intersacc_interval(z)))])

                subplot(2,3,3)
                try
                zp=find(params.saccade.eyedata.arej == 1);
                hist(eyedata.saccdur(zp),0:1:150)
                set(gca,'xlim',[0 150])
                title(['median saccdur: ' num2str(median(eyedata.saccdur(zp))) ', selno: ' num2str(length(zp))])
                catch
                end
                
                subplot(2,3,6)
                try
                hist(eyedata.intersacc_interval(zp),0:50:2000)
                set(gca,'xlim',[0 1000])
                z = find(eyedata.intersacc_interval(zp)<=bool.isi_limit);
                
                title(['median isi ' num2str(median(eyedata.intersacc_interval(zp(z))))])
                catch
                end
                
            end
            
            
            axes('Position',[0 0.98 1 0.2],'Visible','off');
            text(0.5,0,[fname ' - ec_eo_ratio ' num2str(eyedata.ec_eo_ratio), ', number of saccades: ' num2str(size(eyedata.saccades0,1)) '/' num2str(length(zp))],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
            
            aa = sprintf('%d', bool.det_met);
            
            if bool.print == 1
                print ('-djpeg', '-r300', [directory 'eye0' num2str(bool.detection_method) '_' fname '.jpg']);
            end
            
            
            
            if bool.image == 2
                close all
            end
        end
        
    end
end
close (waitb)
