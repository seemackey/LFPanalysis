
function [filenamesout2,filenames_nowave] =  snr_wavelet08cm(directory3,directory4,filenames, wmethod, frqrange,  epoch_tframe_orig, bool,filestart)

selch_ai        = [1:4];
selch_lfp       = [1:4];
a00= 0;

filter_hg       = [70 160];

%%%%%%%%% important: bool.eegc = 1 only saves power, bool.eegc = 2 saves
%%%%%%%%% phase this will be corrected


%%%% examples:

%%%%%%%%% 2020 saccades03 hires-%% Noelle changed it from 5 3 3 to 4 3 3
% close all; clear all; tic
% directory3        = ['D:\dyneye\spont\spont_kieran\spontaneous_addon\spontaneous\x\'];
% directory4        = ['D:\dyneye\spont\spont_kieran\spontaneous_addon\spontaneous\wavelet\'];
% filenames         = {};
% wmethod           = [4 3 3]; %%% 9 7 3 high frq, 4 1 3 low, 5 1 3 all
% frqrange          = [];
% epoch_tframe_orig = [-1000 1000];
% bool.notch60hz    = 1;
% bool.lfp          = 0;
% bool.csd          = 1;
% bool.mua          = 1;
% bool.bipolar      = 0;
% bool.hgamma       = 0;
% bool.arej         = 0;
% bool.eegc          = 0; bool.eegm          = 0; bool.eegb          = 0;
% bool.eege          = 0; bool.eegl          = 0; bool.eega          = 0; bool.eegh          = 0;
% bool.downsample   = 1;
% bool.ph           = 1;
% bool.saccade_arej = 1;
% bool.ttypes       = 0; %%%%%% if 0, all trigger types will be treated as the same,
% %%%%%% if 1, they will be treated as individuals,
% %%%%%% if 2, they will be renumbered based on distance from the deviant
% %%%%%% if 3, epochs will be aligned to pattern onsets
% [filenamesout2,filenames_nowave] =  snr_wavelet08(directory3,directory4,filenames, wmethod, frqrange,  epoch_tframe_orig, bool);
% toc

%%%%%%%%% 2020 saccades03 lowres
% close all; clear all; tic
% directory3        = ['D:\dyneye\spont\spont_kieran\spontaneous_addon\spontaneous\x\'];
% directory4        = ['D:\dyneye\spont\spont_kieran\spontaneous_addon\spontaneous\wavelet\'];
% filenames         = {};
% wmethod           = [2 3 3]; %%% 9 7 3 high frq, 4 1 3 low, 5 1 3 all
% frqrange          = [];
% epoch_tframe_orig = [-1000 1000];
% bool.notch60hz    = 1;
% bool.lfp          = 0;
% bool.csd          = 1;
% bool.mua          = 1;
% bool.bipolar      = 0;
% bool.hgamma       = 0;
% bool.arej         = 0;
% bool.eegc          = 0; bool.eegm          = 0; bool.eegb          = 0;
% bool.eege          = 0; bool.eegl          = 0; bool.eega          = 0; bool.eegh          = 0;
% bool.downsample   = 1;
% bool.ph           = 1;
% bool.saccade_arej = 1;
% bool.ttypes       = 0; %%%%%% if 0, all trigger types will be treated as the same,
% %%%%%% if 1, they will be treated as individuals,
% %%%%%% if 2, they will be renumbered based on distance from the deviant
% %%%%%% if 3, epochs will be aligned to pattern onsets
% [filenamesout2,filenames_nowave] =  snr_wavelet08(directory3,directory4,filenames, wmethod, frqrange,  epoch_tframe_orig, bool);
% toc

%%%%%%%%%% 2020 saccades02
% close all; clear all; tic
% directory3        = ['C:\Data\saccdata2\led_a1_new\contproc\noelle\'];
% directory4        = ['C:\Data\saccdata2\led_a1_new\wavelet\noelle\'];
% filenames         = {};
% wmethod           = [5 3 3]; %%% 9 7 3 high frq, 4 1 3 low, 5 1 3 all
% frqrange          = [0.8 55];
% epoch_tframe_orig = [-500 500];
% bool.notch60hz    = 1;
% bool.lfp          = 0;
% bool.csd          = 1;
% bool.mua          = 1;
% bool.bipolar      = 0;
% bool.hgamma       = 0;
% bool.arej         = 1;
% bool.eegc          = 0; bool.eegm          = 0; bool.eegb          = 0;
% bool.eege          = 0; bool.eegl          = 0; bool.eega          = 0; bool.eegh          = 0;
% bool.downsample   = 1;
% bool.ph           = 1;
% bool.saccade_arej = 1;
% bool.ttypes       = 0; %%%%%% if 0, all trigger types will be treated as the same,
% %%%%%% if 1, they will be treated as individuals,
% %%%%%% if 2, they will be renumbered based on distance from the deviant
% %%%%%% if 3, epochs will be aligned to pattern onsets
% [filenamesout2,filenames_nowave] =  snr_wavelet08(directory3,directory4,filenames, wmethod, frqrange,  epoch_tframe_orig, bool);
% toc

%%%%%%%%%% 2020 saccades01
% close all; clear all; tic
% directory3        = ['C:\Data\saccdata2\led_a1_new\contproc\noelle\'];
% directory4        = ['C:\Data\saccdata2\led_a1_new\wavelet\noelle\'];
% filenames         = {};
% wmethod           = [5 3 3]; %%% 9 7 3 high frq, 4 1 3 low, 5 1 3 all
% frqrange          = [];
% epoch_tframe_orig = [-500 500];
% bool.notch60hz    = 1;
% bool.lfp          = 0;
% bool.csd          = 1;
% bool.mua          = 1;
% bool.bipolar      = 0;
% bool.hgamma       = 0;
% bool.arej         = 1;
% bool.eegc          = 0; bool.eegm          = 0; bool.eegb          = 0;
% bool.eege          = 0; bool.eegl          = 0; bool.eega          = 0; bool.eegh          = 0;
% bool.downsample   = 1;
% bool.ph           = 1;
% bool.saccade_arej = 1;
% bool.ttypes       = 0; %%%%%% if 0, all trigger types will be treated as the same,
% %%%%%% if 1, they will be treated as individuals,
% %%%%%% if 2, they will be renumbered based on distance from the deviant
% %%%%%% if 3, epochs will be aligned to pattern onsets
% [filenamesout2,filenames_nowave] =  snr_wavelet08(directory3,directory4,filenames, wmethod, frqrange,  epoch_tframe_orig, bool);
% toc


downsample_adrate   = 1;

omega               = 6;
% djorig              = 0.0905;   % this was the real original one
djorig              = 0.088;

%%%%%%%%%%%%%%      % 1     % 2     %3      % 4     % 5     % 6     % 7     %8      %9
adrates     = [     50      100     250     500     1000    2500    5000    10000  20000];

%%%%%%%%%%%%%%%%    % 1     % 2     % 3     % 4     % 5     % 6     %7      %8
frqrange1s  = [     0.01    0.5    0.8      1       5       25      65      250];

%%%%%%%%%%%%%%%%% 1         % 2         % 3         % 4           %5            %6
djs         = [ djorig/4  djorig/2     djorig      djorig*2    djorig*4     djorig*8];

%%%%%%%%% to recreate all origianl wavelet files formats, bool.downsample
%%%%%%%%% has to be set to 0, we never downsampled the continuous wavelet
%%%%%%%%% before

%%%%%%%%% @osw0 can be recreated using the settings:  wmethod = [4 3 3];
%%%%%%%%% @osw5 can be recreated using the settings:  wmethod = [2 2 3];

filenames_nowave = {};

if isempty(filenames)
    filenames=[];
    f=dir( [ directory3 '*@os*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
    f=dir( [ directory3 '*@osim.mat']);
    for ii=1:1:size(f,1)
        ff=f(ii).name;
        filenames{i+ii}=ff;
    end
end

[s,mess,messid] = mkdir(directory4);

frqrange_backup = frqrange;

waitb           = waitbar(0);

tic

for filecik=filestart:length(filenames)
    
    wraw=[];
    wlfp=[];
    wai=[];
    
    epoch_tframe = epoch_tframe_orig;
    
    frqrange = frqrange_backup;
    
    wraw.adrate = [];
    wraw.cnte = [];
    wraw.cntm = [];
    wraw.cntb = [];
    wraw.cntc = [];
    wraw.cnth = [];
    wai.adrate = [];
    wai.cnt = [];
    wai.cnt_ph = [];
    wai.cnt_po = [];
    wlfp.adrate = [];
    wlfp.cnt = [];
    wlfp.cnt_ph = [];
    wlfp.cnt_po = [];
    wraw.cntc_po = [];
    wraw.cntc_ph = [];
    wraw.cnte_po = [];
    wraw.cnte_ph = [];
    wraw.cntm_po = [];
    wraw.cntm_ph = [];
    wraw.cntb_po = [];
    wraw.cntb_ph = [];
    wraw.cnth_po = [];
    wraw.cnth_ph = [];
    
    wraw.avg.pob = {};
    wraw.avg.poc = {};
    wraw.avg.pom = {};
    wraw.avg.poe = {};
    wraw.avg.poh = {};
    wlfp.avg.po = {};
    wai.avg.po = {};
    wraw.avg.avg = {};
    wlfp.avg.avg = {};
    wai.avg.avg = {};
    
    
    speaks = [];
    
    waitb           = waitbar(0);
    a2=['loading file ' filenames{filecik}];
    waitbar(0, waitb,a2);
    
    load([directory3 filenames{filecik}])
    
    trig.adrate_orig = craw.adrate;
    
    
    newadrate   = adrates(wmethod(1));
    
    if isempty(frqrange)
        frqrange    = [frqrange1s(wmethod(2)) newadrate*0.45];
    else
        wmethod(2)=length(frqrange1s)+1;
    end
    dj          =  djs(wmethod(3));
    
    
    params.newadrate    = newadrate;
    params.origadrate   = craw.adrate;
    params.frqrange     = frqrange;
    params.dj           = dj;
    params.omega        = omega;
    
    try
        reprate = craw.reprate;
    catch
        reprate = 0;
    end
    
    
    
    %%%%%%%%%%%%%%%% craw downsample
    
    %%%%%% field sua mua and csd
    if ~isempty(craw.adrate) & bool.csd+bool.mua+bool.lfp+bool.bipolar+bool.hgamma>0
        
        wraw.reprate = reprate;
        
        a2=['transforming raw data ' filenames{filecik}];
        waitbar(0, waitb,a2);
        
        
        
        
        % 60 Hz filter
        
        if bool.notch60hz == 1
            d = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                'DesignMethod','butter','SampleRate',craw.adrate);
            for i=1:size(craw.cnt,1)
                craw.cnt(i,:)=filtfilt(d,craw.cnt(i,:));
            end
            
            %%%%%%%%%% old method
            %             wo = 60/(adrate/2);
            %             bw = wo/100;
            %             [b,a] = iirnotch(wo,bw);
            %             for i=1:size(cnt,1)
            %                 cnt(i,:)=filtfilt(b,a,cnt(i,:));
            %             end
            
            if ~isempty(clfp.adrate)
                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                    'DesignMethod','butter','SampleRate',clfp.adrate);
                for i=1:size(clfp.cnt,1)
                    clfp.cnt(i,:)=filtfilt(d,clfp.cnt(i,:));
                end
            end
            
        end
        
        adrate = craw.adrate;
        
        
        
        wraw.adrate_orig = adrate;
        wraw.adrate = newadrate;
        wraw.wmethod = wmethod;
        wraw.frqrange = frqrange;
        
        %%%% new filter to get rid of "DC shift"
        %         n = 2;
        %         Wn = [0.2 9000]/(adrate/2);
        %         [b,a] = butter(n,Wn);
        %         for i=1:size(cnt,1)
        %             cnt(i,:)=filtfilt(b,a,cnt(i,:));
        %         end
        
        %                  for i=1:size(cnt,1)
        %              cnt(i,:)=cnt(i,:)-mean(cnt(i,:));
        %          end
        
        
        
        if frqrange(2)*1.25>300
            filtere = [frqrange(1)/2 newadrate*0.49];
        else
            filtere = [frqrange(1)/2 300];
        end
        filteru = [300 6000];
        filtertype = 1;
        
        
        params.filtere = filtere;
        params.filteru = filteru;
        params.filter_hg = filter_hg;
        params.dj = dj;
        params.omega = omega;
        
        [~, cnte, cntm, cntc, ~, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);
        %%%%%% create high gamma cnt signal
        if bool.hgamma          == 1
            cnth = zeros(size(cntb));
            n = 2;
            Wn = filter_hg/(newadrate/2);
            [b,a] = butter(n,Wn);
            for i1=1:size(cntb,1)
                cnth(i1,:)=filtfilt(b,a,cntb(i1,:));
            end
            for i1=1:size(cnth,1)
                cnth(i1,:)=abs(hilbert(cnth(i1,:)));
            end
        end
        
        %         [~, cnte, cntm, cntc, ~, cntb, cnth] = module_cnt05b(craw, newadrate, filtere, filteru, filter_hg, filtertype);
        %         % rectifying hgamma
        %         for i1=1:size(cnth,1)
        %             cnth(i1,:)=abs(hilbert(cnth(i1,:)));
        %         end
        
        
        
        cntm(end+1,:)=mean(cntm,1);
        cnte(end+1,:)=mean(abs(hilbert(cnte)),1);
        cntc(end+1,:)=mean(abs(hilbert(cntc)),1);
        cntb(end+1,:)=mean(abs(hilbert(cntb)),1);
        if bool.hgamma          == 1
            cnth(end+1,:)=mean(abs(hilbert(cnth)),1);
        end
        
        adrate = newadrate;
        
        
        
        if bool.lfp == 1
            wraw.cnte = cnte;
        end
        if bool.csd == 1
            wraw.cntc = cntc;
        end
        if bool.mua == 1
            wraw.cntm = cntm;
        end
        if bool.bipolar == 1
            wraw.cntb = cntb;
        end
        if bool.hgamma == 1
            wraw.cnth = cnth;
        end
        
        clear cnte cntc cntm cntb cnth
        
        
        
        %%%% wavelet
        
        
        if bool.lfp == 1
            for chancik=1:size(wraw.cnte,1)
                a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', lfp'];
                waitbar(chancik/size(wraw.cnte,1), waitb,a2);
                [wave,period,scale,coi] = basewave6(wraw.cnte(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
                
                
                
                if chancik == 1
                    
                    wraw.cnte_po = zeros(size(wraw.cnte,1),length(period),size(wraw.cnte,2));
                end
                
                pox=(((abs(wave)/24).^2)*1000)/adrate;
                for i=1:size(pox,1)
                    pox(i,:)=sqrt(pox(i,:)./period(i));
                end
                wraw.cnte_po(chancik,:,:)=pox;
                
                if bool.ph == 1
                    if chancik == 1
                        wraw.cnte_ph = zeros(size(wraw.cnte,1),length(period),size(wraw.cnte,2));
                        
                    end
                    wraw.cnte_ph(chancik,:,:)=angle(wave);
                    
                    
                end
                
            end
        end
        
        if bool.csd == 1
            for chancik=1:size(wraw.cntc,1)
                a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', csd'];
                waitbar(chancik/size(wraw.cntc,1), waitb,a2);
                [wave,period,scale,coi] = basewave6(wraw.cntc(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
                
                
                if chancik == 1
                    
                    wraw.cntc_po = zeros(size(wraw.cntc,1),length(period),size(wraw.cntc,2));
                end
                
                pox=(((abs(wave)/24).^2)*1000)/adrate;
                for i=1:size(pox,1)
                    pox(i,:)=sqrt(pox(i,:)./period(i));
                end
                wraw.cntc_po(chancik,:,:)=pox;
                
                if bool.ph == 1
                    if chancik == 1
                        wraw.cntc_ph = zeros(size(wraw.cntc,1),length(period),size(wraw.cntc,2));
                        
                    end
                    wraw.cntc_ph(chancik,:,:)=angle(wave);
                end
                
            end
            
        end
        
        if bool.mua == 1
            for chancik=1:size(wraw.cntm,1)
                a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', mua'];
                waitbar(chancik/size(wraw.cntm,1), waitb,a2);
                [wave,period,scale,coi] = basewave6(wraw.cntm(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
                if chancik == 1
                    
                    wraw.cntm_po = zeros(size(wraw.cntm,1),length(period),size(wraw.cntm,2));
                end
                
                pox=(((abs(wave)/24).^2)*1000)/adrate;
                for i=1:size(pox,1)
                    pox(i,:)=sqrt(pox(i,:)./period(i));
                end
                wraw.cntm_po(chancik,:,:)=pox;
                
                if bool.ph == 1
                    if chancik == 1
                        wraw.cntm_ph = zeros(size(wraw.cntm,1),length(period),size(wraw.cntm,2));
                        
                    end
                    wraw.cntm_ph(chancik,:,:)=angle(wave);
                end
                
            end
        end
        
        if bool.bipolar == 1
            for chancik=1:size(wraw.cntb,1)
                a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', bipolar'];
                waitbar(chancik/size(wraw.cntb,1), waitb,a2);
                [wave,period,scale,coi] = basewave6(wraw.cntb(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
                if chancik == 1
                    
                    wraw.cntb_po = zeros(size(wraw.cntb,1),length(period),size(wraw.cntb,2));
                end
                
                pox=(((abs(wave)/24).^2)*1000)/adrate;
                for i=1:size(pox,1)
                    pox(i,:)=sqrt(pox(i,:)./period(i));
                end
                wraw.cntb_po(chancik,:,:)=pox;
                
                if bool.ph == 1
                    if chancik == 1
                        wraw.cntb_ph = zeros(size(wraw.cntb,1),length(period),size(wraw.cntb,2));
                        
                    end
                    wraw.cntb_ph(chancik,:,:)=angle(wave);
                end
                
            end
        end
        
        if bool.hgamma == 1
            
            for chancik=1:size(wraw.cnth,1)
                a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', hgamma'];
                waitbar(chancik/size(wraw.cnth,1), waitb,a2);
                [wave,period,scale,coi] = basewave6(wraw.cnth(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
                if chancik == 1
                    
                    wraw.cnth_po = zeros(size(wraw.cnth,1),length(period),size(wraw.cnth,2));
                end
                
                pox=(((abs(wave)/24).^2)*1000)/adrate;
                for i=1:size(pox,1)
                    pox(i,:)=sqrt(pox(i,:)./period(i));
                end
                wraw.cnth_po(chancik,:,:)=pox;
                
                if bool.ph == 1
                    if chancik == 1
                        wraw.cnth_ph = zeros(size(wraw.cnth,1),length(period),size(wraw.cnth,2));
                        
                    end
                    wraw.cnth_ph(chancik,:,:)=angle(wave);
                end
                
            end
        end
        
        wraw.frq=1./period;
        wraw.frqrange = frqrange;
        
        
        
        dt              = 1000/wraw.adrate;
        wraw.time      = (0:size(wraw.cntc,2)-1).*dt;
        
        clear wave
        
    end
    
    
    
    %%%%%%%%% trigger downsample
    
    for i1=1:length(trig.anatrig)
        if ~isempty(trig.anatrig{i1})
            trig.anatrig{i1}    = round(trig.anatrig{i1}./(trig.adrate/newadrate));
        end
    end
    for i1=1:size(trig.digtrig,2)
        trig.digtrig(:,i1)    = round(trig.digtrig(:,i1)./(trig.adrate/newadrate));
    end
    try
        trig.adrate_orig = trig.adrate;
        trig.adrate = newadrate;
    catch
    end
    
    %%%%%%%%%%%%%%%% cai down or upsample
    if ~isempty(cai.adrate) && ~isempty(selch_ai)
        
        
        wai.reprate = reprate;
        
        
        
        z=find(selch_ai)>size(cai.cnt,1);
        selch_ai(z)=[];
        
        
        adrate = cai.adrate;
        xx = cai.cnt(selch_ai,:);
        clear cai
        
        if newadrate>adrate
            xx          = resample(xx',newadrate,adrate)';
        elseif newadrate==adrate
            
        else
            
            % filtering before downsample
            
            nyq = newadrate/2;
            n = 2;
            Wn = nyq/(adrate/2);
            [b,a] = butter(n,Wn,'low');
            for i=1:size(xx,1)
                xx(i,:)=filtfilt(b,a,xx(i,:));
            end
            
            % downsample or resmaple
            downsampleby    = adrate/newadrate;
            if downsampleby-round(downsampleby)==0
                xx          = downsample(xx',downsampleby)';
            else
                xx          = resample(xx',newadrate,adrate)';
            end
            
        end
        
        
        wai.cnt         = xx;
        wai.selch       = selch_ai;
        wai.adrate      = newadrate;
        wai.adrate_orig = adrate;
        
        cnt     = wai.cnt;
        adrate  = wai.adrate;
        
        
        for chancik=1:size(cnt,1)
            
            a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', ai'];
            waitbar(chancik/size(cnt,1), waitb,a2);
            
            [wave,period,scale,coi] = basewave6(cnt(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
            
            if chancik == 1
                wai.cnt_po = zeros(size(wai.cnt,1),length(period),size(wai.cnt,2));
            end
            pox=(((abs(wave)/24).^2)*1000)/adrate;
            for i=1:size(pox,1)
                pox(i,:)=sqrt(pox(i,:)./period(i));
            end
            wai.cnt_po(chancik,:,:)=pox;
            
            if bool.ph == 1
                if chancik == 1
                    wai.cnt_ph = zeros(size(wai.cnt,1),length(period),size(wai.cnt,2));
                end
                wai.cnt_ph(chancik,:,:)=angle(wave);
            end
            
        end
        wai.frq=1./period;
        wai.frqrange = frqrange;
        wai.wmethod = wmethod;
        
        dt              = 1000/wai.adrate;
        wai.time        = (0:size(wai.cnt,2)-1).*dt;
        
    end
    
    %%%%%%%%%%%%%%%% clfp down or upsample
    if ~isempty(clfp.adrate) && ~isempty(selch_lfp)
        
        wlfp.reprate = reprate;
        
        z=find(selch_lfp)>size(clfp.cnt,1);
        selch_lfp(z)=[];
        
        adrate = clfp.adrate;
        xx = clfp.cnt(selch_lfp,:);
        clear clfp
        
        if newadrate>adrate
            xx          = resample(xx',newadrate,adrate)';
        else
            
            % filtering before downsample
            nyq = newadrate/2;
            n = 2;
            Wn = nyq/(adrate/2);
            [b,a] = butter(n,Wn,'low');
            for i=1:size(xx,1)
                xx(i,:)=filtfilt(b,a,xx(i,:));
            end
            
            % downsample or resmaple
            downsampleby    = adrate/newadrate;
            if downsampleby-round(downsampleby)==0
                xx          = downsample(xx',downsampleby)';
            else
                xx          = resample(xx',newadrate,adrate)';
            end
            
        end
        
        wlfp.cnt         = xx;
        wlfp.selch       = selch_lfp;
        wlfp.adrate      = newadrate;
        wlfp.adrate_orig = adrate;
        
        cnt     = wlfp.cnt;
        adrate  = wlfp.adrate;
        
        
        
        
        
        
        
        for chancik=1:size(cnt,1)
            
            a2=[num2str(filecik) ' - ' filenames{filecik} ',ch ' num2str(chancik) ', lfp'];
            waitbar(chancik/size(cnt,1), waitb,a2);
            
            [wave,period,scale,coi] = basewave6(cnt(chancik,:),adrate,frqrange(1),frqrange(2),omega,dj,0,1);
            
            
            if chancik == 1
                wlfp.cnt_po = zeros(size(wlfp.cnt,1),length(period),size(wlfp.cnt,2));
            end
            pox=(((abs(wave)/24).^2)*1000)/adrate;
            for i=1:size(pox,1)
                pox(i,:)=sqrt(pox(i,:)./period(i));
            end
            wlfp.cnt_po(chancik,:,:)=pox;
            if bool.ph == 1
                if chancik == 1
                    wlfp.cnt_ph = zeros(size(wlfp.cnt,1),length(period),size(wlfp.cnt,2));
                end
                wlfp.cnt_ph(chancik,:,:)=angle(wave);
            end
        end
        wlfp.frq=1./period;
        wlfp.frqrange = frqrange;
        wlfp.wmethod = wmethod;
        dt              = 1000/wlfp.adrate;
        wlfp.time        = (0:size(wlfp.cnt,2)-1).*dt;
    end
    
    
    
    clear wave pox
    % ylabels     = num2str(frq','%11.3g');
    
    a = '';
    b = '';
    
    if bool.downsample == 1
        a = [a 'd'];
    end
    if bool.eegc > 0 | bool.eegl > 0 | bool.eega > 0 | bool.eegm > 0 | bool.eegb > 0 | bool.eege > 0
        a = [a 'e'];
    end
    
    
    
    wmethod_number = wmethod(1)*100+wmethod(2)*10+wmethod(3);
    
    findex2 = [num2str(wmethod_number) a b];
    
    if  bool.saccade_arej == 1
        aa = 'sa';
    else
        aa = '';
    end

    fn = filenames{filecik};
    k = strfind(fn,'@os');
    fn_new = [fn(1:k-1) '@osw' findex2 aa fn(k+3:end)];
    filenamesout2{filecik}=fn_new;
    
    fn = filenamesout2{filecik};

    [wraw,wlfp,wai, epoch_tframe] = module_wavelet_avg06(wraw, wlfp, wai, trig, params, epoch_tframe, bool);
    wraw.avg.adrate = wraw.adrate;
    wai.avg.adrate = wai.adrate;
    wlfp.avg.adrate = wlfp.adrate;
    
    
    close (waitb)
    
    if size(wraw.cntc_po,2)==0
        a00= a00+1;
        filenames_nowave(a00) = filenames(filecik);
    else
        
        
        
        %[speaks.peaks,speaks.surfcohere,speaks.ccoef,speaks.params,speaks.posall] = module_spectral_peaks02(wraw, wlfp);
        
        if bool.downsample == 1
            wraw.cntc_ph = [];
            wraw.cntm_ph = [];
            wraw.cntb_ph = [];
            wraw.cnte_ph = [];
            wraw.cnth_ph = [];
            wai.cnt_ph = [];
            wlfp.cnt_ph = [];
            newadrate = downsample_adrate;
            nyq = newadrate/2.01;
            
            adrate  = wraw.adrate;
            if isempty(adrate)
                adrate  = wlfp.adrate;
            end
            if isempty(adrate)
                adrate  = wai.adrate;
            end
            
            %%%% lowpass filter
            n = 2;
            Wn = nyq/(adrate/2);
            [b,a] = butter(n,Wn,'low');
            for i1=1:size(wraw.cntc,1)
                wraw.cntc(i1,:)=filtfilt(b,a,wraw.cntc(i1,:));
            end
            for i1=1:size(wraw.cntm,1)
                wraw.cntm(i1,:)=filtfilt(b,a,wraw.cntm(i1,:));
            end
            for i1=1:size(wraw.cntb,1)
                wraw.cntb(i1,:)=filtfilt(b,a,wraw.cntb(i1,:));
            end
            for i1=1:size(wraw.cnte,1)
                wraw.cnte(i1,:)=filtfilt(b,a,wraw.cnte(i1,:));
            end
            for i1=1:size(wraw.cnth,1)
                wraw.cnth(i1,:)=filtfilt(b,a,wraw.cnth(i1,:));
            end
            
            for i1=1:size(wai.cnt,1)
                wai.cnt(i1,:)=filtfilt(b,a,wai.cnt(i1,:));
            end
            for i1=1:size(wlfp.cnt,1)
                wlfp.cnt(i1,:)=filtfilt(b,a,wlfp.cnt(i1,:));
            end
            
            
            for i1=1:size(wraw.cntc_po,1)
                for i2=1:size(wraw.cntc_po,2)
                    wraw.cntc_po(i1,i2,:)=filtfilt(b,a,squeeze(wraw.cntc_po(i1,i2,:)));
                end
            end
            for i1=1:size(wraw.cntm_po,1)
                for i2=1:size(wraw.cntm_po,2)
                    wraw.cntm_po(i1,i2,:)=filtfilt(b,a,squeeze(wraw.cntm_po(i1,i2,:)));
                end
            end
            for i1=1:size(wraw.cntb_po,1)
                for i2=1:size(wraw.cntb_po,2)
                    wraw.cntb_po(i1,i2,:)=filtfilt(b,a,squeeze(wraw.cntb_po(i1,i2,:)));
                end
            end
            for i1=1:size(wraw.cnte_po,1)
                for i2=1:size(wraw.cnte_po,2)
                    wraw.cnte_po(i1,i2,:)=filtfilt(b,a,squeeze(wraw.cnte_po(i1,i2,:)));
                end
            end
            for i1=1:size(wraw.cnth_po,1)
                for i2=1:size(wraw.cnth_po,2)
                    wraw.cnth_po(i1,i2,:)=filtfilt(b,a,squeeze(wraw.cnth_po(i1,i2,:)));
                end
            end
            
            for i1=1:size(wai.cnt_po,1)
                for i2=1:size(wai.cnt_po,2)
                    wai.cnt_po(i1,i2,:)=filtfilt(b,a,squeeze(wai.cnt_po(i1,i2,:)));
                end
            end
            for i1=1:size(wlfp.cnt_po,1)
                for i2=1:size(wlfp.cnt_po,2)
                    wlfp.cnt_po(i1,i2,:)=filtfilt(b,a,squeeze(wlfp.cnt_po(i1,i2,:)));
                end
            end
            
            
            
            
            % downsample or resmaple
            downsampleby    = adrate/newadrate;
            
            if bool.csd == 1
                wraw.cntc           = downsample(wraw.cntc',downsampleby)';
            end
            if bool.mua == 1
                wraw.cntm           = downsample(wraw.cntm',downsampleby)';
            end
            if bool.bipolar == 1
                wraw.cntb           = downsample(wraw.cntb',downsampleby)';
            end
            if bool.lfp == 1
                wraw.cnte           = downsample(wraw.cnte',downsampleby)';
            end
            if bool.hgamma == 1
                wraw.cnth           = downsample(wraw.cnth',downsampleby)';
            end
            if ~isempty(wai.cnt)
                wai.cnt             = downsample(wai.cnt',downsampleby)';
            end
            if ~isempty(wlfp.cnt)
                wlfp.cnt            = downsample(wlfp.cnt',downsampleby)';
            end
            
            
            if ~isempty(wai.adrate)
                
                
                
                s3 = length(downsample(squeeze(wai.cnt_po(1,1,:)),downsampleby));
                cntai_po2 = zeros(size(wai.cnt_po,1),size(wai.cnt_po,2),s3);
                
                
                for i1=1:size(wai.cnt_po,1)
                    for i2=1:size(wai.cnt_po,2)
                        cntai_po2(i1,i2,:)=downsample(squeeze(wai.cnt_po(i1,i2,:)),downsampleby);
                    end
                end
                wai.cnt_po      = cntai_po2;
                clear cntai_po2
                wai.adrate = newadrate;
                wai.time = downsample(wai.time,downsampleby);
            end
            if ~isempty(wlfp.adrate)
                s3 = length(downsample(squeeze(wlfp.cnt_po(1,1,:)),downsampleby));
                cntlfp_po2 = zeros(size(wlfp.cnt_po,1),size(wlfp.cnt_po,2),s3);
                for i1=1:size(wlfp.cnt_po,1)
                    for i2=1:size(wlfp.cnt_po,2)
                        cntlfp_po2(i1,i2,:)=downsample(squeeze(wlfp.cnt_po(i1,i2,:)),downsampleby);
                    end
                end
                wlfp.cnt_po     = cntlfp_po2;
                clear cntlfp_po2
                wlfp.adrate = newadrate;
                wlfp.time = downsample(wlfp.time,downsampleby);
            end
            
            if bool.csd == 1
                s3 = length(downsample(squeeze(wraw.cntc_po(1,1,:)),downsampleby));
                cntc_po2 = zeros(size(wraw.cntc_po,1),size(wraw.cntc_po,2),s3);
                for i1=1:size(wraw.cntc_po,1)
                    for i2=1:size(wraw.cntc_po,2)
                        cntc_po2(i1,i2,:)=downsample(squeeze(wraw.cntc_po(i1,i2,:)),downsampleby);
                    end
                end
                wraw.cntc_po    = cntc_po2;
                clear cntc_po2
                
            end
            
            if bool.mua == 1
                s3 = length(downsample(squeeze(wraw.cntm_po(1,1,:)),downsampleby));
                cntm_po2 = zeros(size(wraw.cntm_po,1),size(wraw.cntm_po,2),s3);
                for i1=1:size(wraw.cntm_po,1)
                    for i2=1:size(wraw.cntm_po,2)
                        cntm_po2(i1,i2,:)=downsample(squeeze(wraw.cntm_po(i1,i2,:)),downsampleby);
                    end
                end
                wraw.cntm_po    = cntm_po2;
                clear cntm_po2
                
            end
            
            if bool.bipolar == 1
                s3 = length(downsample(squeeze(wraw.cntb_po(1,1,:)),downsampleby));
                cntb_po2 = zeros(size(wraw.cntb_po,1),size(wraw.cntb_po,2),s3);
                for i1=1:size(wraw.cntb_po,1)
                    for i2=1:size(wraw.cntb_po,2)
                        cntb_po2(i1,i2,:)=downsample(squeeze(wraw.cntb_po(i1,i2,:)),downsampleby);
                    end
                end
                wraw.cntb_po    = cntb_po2;
                clear cntb_po2
                
            end
            
            
            if bool.lfp == 1
                s3 = length(downsample(squeeze(wraw.cnte_po(1,1,:)),downsampleby));
                cnte_po2 = zeros(size(wraw.cnte_po,1),size(wraw.cnte_po,2),s3);
                for i1=1:size(wraw.cnte_po,1)
                    for i2=1:size(wraw.cnte_po,2)
                        cnte_po2(i1,i2,:)=downsample(squeeze(wraw.cnte_po(i1,i2,:)),downsampleby);
                    end
                end
                wraw.cnte_po    = cnte_po2;
                clear cnte_po2
                
            end
            
            if bool.hgamma == 1
                s3 = length(downsample(squeeze(wraw.cnth_po(1,1,:)),downsampleby));
                cnth_po2 = zeros(size(wraw.cnth_po,1),size(wraw.cnth_po,2),s3);
                for i1=1:size(wraw.cnth_po,1)
                    for i2=1:size(wraw.cnth_po,2)
                        cnth_po2(i1,i2,:)=downsample(squeeze(wraw.cnth_po(i1,i2,:)),downsampleby);
                    end
                end
                wraw.cnth_po    = cnth_po2;
                clear cnth_po2
                
            end
            
            if  ~isempty(wraw.adrate)
                wraw.time = downsample(wraw.time,downsampleby);
                wraw.adrate = newadrate;
            end
            
            
            
        else
            
            
            
            
            if ~isempty(wai.adrate)
                cntai_po2 = [];
                for i1=1:size(wai.cnt_po,1)
                    for i2=1:size(wai.cnt_po,2)
                        cntai_po2(i1,i2,:)=resample(squeeze(wai.cnt_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wai.cnt_po      = cntai_po2;
                clear cntai_po2
                wai.adrate = newadrate;
                wai.cnt            = resample(wai.cnt',newadrate,adrate)';
            end
            if ~isempty(wlfp.adrate)
                cntlfp_po2 = [];
                for i1=1:size(wlfp.cnt_po,1)
                    for i2=1:size(wlfp.cnt_po,2)
                        cntlfp_po2(i1,i2,:)=resample(squeeze(wlfp.cnt_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wlfp.cnt_po     = cntlfp_po2;
                clear cntlfp_po2
                wlfp.adrate = newadrate;
                wlfp.cnt            = resample(wlfp.cnt',newadrate,adrate)';
            end
            
            
            if bool.csd == 1
                cntc_po2 = [];
                for i1=1:size(wraw.cntc_po,1)
                    for i2=1:size(wraw.cntc_po,2)
                        cntc_po2(i1,i2,:)=resample(squeeze(wraw.cntc_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wraw.cntc_po    = cntc_po2;
                clear cntc_po2
                wraw.cntc            = resample(wraw.cntc',newadrate,adrate)';
            end
            
            if bool.mua == 1
                cntm_po2 = [];
                for i1=1:size(wraw.cntm_po,1)
                    for i2=1:size(wraw.cntm_po,2)
                        cntm_po2(i1,i2,:)=resample(squeeze(wraw.cntm_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wraw.cntm_po    = cntm_po2;
                clear cntm_po2
            end
            
            if bool.bipolar == 1
                cntb_po2 = [];
                for i1=1:size(wraw.cntb_po,1)
                    for i2=1:size(wraw.cntb_po,2)
                        cntb_po2(i1,i2,:)=resample(squeeze(wraw.cntb_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wraw.cntb_po    = cntb_po2;
                clear cntb_po2
            end
            
            if bool.lfp == 1
                cnte_po2 = [];
                for i1=1:size(wraw.cnte_po,1)
                    for i2=1:size(wraw.cnte_po,2)
                        cnte_po2(i1,i2,:)=resample(squeeze(wraw.cnte_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wraw.cnte_po    = cnte_po2;
                clear cnte_po2
            end
            
            if bool.hgamma == 1
                cnth_po2 = [];
                for i1=1:size(wraw.cnth_po,1)
                    for i2=1:size(wraw.cnth_po,2)
                        cnth_po2(i1,i2,:)=resample(squeeze(wraw.cnth_po(i1,i2,:)),newadrate,adrate);
                    end
                end
                wraw.cnth_po    = cnth_po2;
                clear cnth_po2
            end
            
            wraw.time = resample(wraw.time,newadrate,adrate);
            wraw.adrate = newadrate;
            
            
        end
        
        %%%%%%%%% trigger downsample
        for i1=1:length(trig.anatrig)
            if ~isempty(trig.anatrig{i1})
                trig.anatrig{i1}    = round(trig.anatrig{i1}./(trig.adrate/newadrate));
            end
        end
        for i1=1:size(trig.digtrig,2)
            trig.digtrig(:,i1)    = round(trig.digtrig(:,i1)./(trig.adrate/newadrate));
        end
        try
            trig.adrate = newadrate;
        catch
        end
    end
    
    params.bool = bool;
    
    save([directory4 filenamesout2{filecik}], 'wraw','wai','wlfp','trig','otherdata','params','-mat','-v7.3') % removed save of speaks 2024, chase
    
    
    
    %%%%%%%%%
    %%%%%%%%% images
    
    if bool.images == 1
        
        mapvars_po2_eeg ={};
        mapvars_itc2_eeg ={};
        
        axposes1{1} = [0.05 0.52 0.94 0.42;  0.05 0.05 0.94 0.4];
        axposes1{2} = [0.05 0.52 0.44 0.44;  0.05 0.27 0.94 0.18; 0.55 0.52 0.44 0.44; 0.05 0.05 0.94 0.18];
        axposes1{3} = [0.05 0.62 0.28 0.34;  0.05 0.41 0.94 0.14; 0.38 0.62 0.28 0.34;  0.05 0.23 0.94 0.14; 0.71 0.62 0.28 0.34; 0.05 0.05 0.94 0.14];
        axposes1{4} = [0.03 0.62 0.21 0.32;  0.05 0.47 0.94 0.105; 0.28 0.62 0.21 0.32;  0.05 0.33 0.94 0.105; 0.53 0.62 0.21 0.32; 0.05 0.19 0.94 0.105; 0.78 0.62 0.21 0.32; 0.05 0.05 0.94 0.105];
        axposes1{5} = [0.03 0.62 0.21 0.32;  0.05 0.47 0.94 0.105; 0.28 0.62 0.21 0.32;  0.05 0.33 0.94 0.105; 0.53 0.62 0.21 0.32; 0.05 0.19 0.94 0.105; 0.78 0.62 0.21 0.32; 0.05 0.05 0.94 0.105];
        
        
        % first get all power (po) and itc values into cells 
        if ~isempty(wraw.adrate)
            
            imagecount = 0;
            imagecount2 = 0;
            if ~isempty(wraw.cntc_po)
                imagecount=imagecount+1;
                mapvars_po1{imagecount}=squeeze(mean(wraw.cntc_po(1:end-1,:,:),3));
                mapvars_po2{imagecount}=squeeze(mean(wraw.cntc_po(1:end-1,:,:),1));
                if ~isempty(wraw.avg.poc)
                    imagecount2=imagecount2+1;
                    for trigcik = 1:length(wraw.avg.poc)
                        mapvars_po1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poc{trigcik}(1,1:end-1,:,:),4));
                        mapvars_itc1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcc{trigcik}(1,1:end-1,:,:),4));
                        mapvars_po2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poc{trigcik}(1,1:end-1,:,:),2));
                        mapvars_itc2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcc{trigcik}(1,1:end-1,:,:),2));
                    end
                    
                end
                ilabels{imagecount}  = 'CSD';
            end
            if ~isempty(wraw.cntm_po)
                imagecount=imagecount+1;
                mapvars_po1{imagecount}=squeeze(mean(wraw.cntm_po(1:end-1,:,:),3));
                mapvars_po2{imagecount}=squeeze(mean(wraw.cntm_po(1:end-1,:,:),1));
                if ~isempty(wraw.avg.pom)
                    imagecount2=imagecount2+1;
                    for trigcik = 1:length(wraw.avg.pom)
                        mapvars_po1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.pom{trigcik}(1,1:end-1,:,:),4));
                        mapvars_itc1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcm{trigcik}(1,1:end-1,:,:),4));
                        mapvars_po2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.pom{trigcik}(1,1:end-1,:,:),2));
                        mapvars_itc2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcm{trigcik}(1,1:end-1,:,:),2));
                    end
                end
                ilabels{imagecount}  = 'MUA';
            end
            if ~isempty(wraw.cntb_po)
                imagecount=imagecount+1;
                mapvars_po1{imagecount}=squeeze(mean(wraw.cntb_po(1:end-1,:,:),3));
                mapvars_po2{imagecount}=squeeze(mean(wraw.cntb_po(1:end-1,:,:),1));
                if ~isempty(wraw.avg.pob)
                    imagecount2=imagecount2+1;
                    for trigcik = 1:length(wraw.avg.pob)
                        mapvars_po1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.pob{trigcik}(1,1:end-1,:,:),4));
                        mapvars_itc1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcb{trigcik}(1,1:end-1,:,:),4));
                        mapvars_po2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.pob{trigcik}(1,1:end-1,:,:),2));
                        mapvars_itc2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itcb{trigcik}(1,1:end-1,:,:),2));
                    end
                end
                ilabels{imagecount}  = 'Bip';
            end
            if ~isempty(wraw.cnte_po)
                imagecount=imagecount+1;
                mapvars_po1{imagecount}=squeeze(mean(wraw.cnte_po(1:end-1,:,:),3));
                mapvars_po2{imagecount}=squeeze(mean(wraw.cnte_po(1:end-1,:,:),1));
                if ~isempty(wraw.avg.poe)
                    imagecount2=imagecount2+1;
                    for trigcik = 1:length(wraw.avg.poe)
                        mapvars_po1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poe{trigcik}(1,1:end-1,:,:),4));
                        mapvars_itc1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itce{trigcik}(1,1:end-1,:,:),4));
                        mapvars_po2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poe{trigcik}(1,1:end-1,:,:),2));
                        mapvars_itc2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itce{trigcik}(1,1:end-1,:,:),2));
                    end
                end
                ilabels{imagecount}  = 'LFP';
            end
            if ~isempty(wraw.cnth_po)
                %%%%%%%%%%%%%%% fix this, temproray
                %%%%%% imagecount=imagecount+1;
                mapvars_po1{imagecount}=squeeze(mean(wraw.cnth_po(1:end-1,:,:),3));
                mapvars_po2{imagecount}=squeeze(mean(wraw.cnth_po(1:end-1,:,:),1));
                if ~isempty(wraw.avg.poh)
                    imagecount2=imagecount2+1;
                    for trigcik = 1:length(wraw.avg.poh)
                        mapvars_po1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poh{trigcik}(1,1:end-1,:,:),4));
                        mapvars_itc1_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itch{trigcik}(1,1:end-1,:,:),4));
                        mapvars_po2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.poh{trigcik}(1,1:end-1,:,:),2));
                        mapvars_itc2_eeg{trigcik,imagecount2}=squeeze(mean(wraw.avg.itch{trigcik}(1,1:end-1,:,:),2));
                    end
                end
                ilabels{imagecount}  = 'HGamma';
            end
            
            ylabels = num2str(wraw.frq','%11.3g');
            dt      = 1000/wraw.adrate;
            time    = wraw.time;
            
            
            
            
            curfig = figure;
            set(curfig,'position',[100   50   1200   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
            colormap jet
            ylabelres=10;
            fsize = 5;
            
            
            
            ic1 = 0;
            
            
            % loop through diff measures (lfp, csd, mua etc.) and plot
            % power and itc
            for icik = 1:imagecount
                
                
                
                
                ic1=ic1+1;
                axpos=axposes1{imagecount}(ic1,:);
                figureax1a=axes('Position',axpos);
                

                mapvar = mapvars_po1{icik};
                surface(1:size(mapvar,2),1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
                set(gca,'xlim',[1 size(mapvar,2)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out','ydir','reverse')
                set(gca,'xlim',[1 length(wraw.frq)],'xtick',1:ylabelres:size(mapvar,2),'xticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                
                cax1a=caxis;
                title([ 'laminar ' ilabels{icik} ' spectrogram,   cax: ' num2str(cax1a(1)) ' - ' num2str(cax1a(2))],'FontSize',8)
                
                ic1=ic1+1;
                axpos=axposes1{imagecount}(ic1,:);
                figureax1a=axes('Position',axpos);
                
                mapvar = mapvars_po2{icik};
                surface(time,1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
                set(gca,'xlim',[time(1) time(end)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out')
                if icik<imagecount
                    set(gca,'xticklabel','')
                end
                set(gca,'ylim',[1 length(wraw.frq)],'ytick',1:ylabelres:size(mapvar,1),'yticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                caxis(cax1a)
                cax3=caxis;
                title([ 'time-frequency ' ilabels{icik} ',   cax: ' num2str(cax3(1)) ' - ' num2str(cax3(2))],'FontSize',8)
            end
            
            % save fig and print fig UNDER CONSTRUCTION %%%%%%%%%%%%%%%%
            fname = filenamesout2{filecik}(1:end-4);
            axes('Position',[0 0.98 1 0.2],'Visible','off');
            text(0.5,0,[fname ''],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
            
            print ('-djpeg', '-r1200', [directory4 'w01_' filenamesout2{filecik}(1:end-4) '.jpg']);
            
            savefig(fname)
            close all
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if imagecount2>0
                
                %%%%%%%%%%% temp fix
                if imagecount2>4
                    imagecount2 = 4;
                end
                
                curfig = figure;
                set(curfig,'position',[100   50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                
                t3 = wraw.avg.time;
                frq = wraw.frq;
                ylabels     = num2str(frq','%11.3g');
                colormap jet
                ylabelres=10;
                fsize = 5;
                
                for trigcik = 1:size(mapvars_po2_eeg,1)
                    ic1 = 0;
                    for icik = 1:imagecount2
                        
                        ic1 = ic1+1;
                        subplot(imagecount2,4,ic1)
                        mapvar=mapvars_po2_eeg{trigcik,icik};
                        surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                        set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                        set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                        caxset1=caxis;
                        caxset1(2)=max(max(mapvar))*0.9;
                        try
                        caxis(caxset1)
                        catch
                        end
                        title (['Amp ' ilabels{icik} ', ' num2str(wraw.avg.sweepno{trigcik}(1)) ' epochs, scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                        
                        ic1 = ic1+1;
                        subplot(imagecount2,4,ic1)
                        mapvar=mapvars_po1_eeg{trigcik,icik};
                        surface(1:size(mapvar,2),1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                        set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                        set(gca,'xlim',[1 size(mapvar,2)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out','ydir','reverse')
                        set(gca,'xlim',[1 length(wraw.frq)],'xtick',1:ylabelres:size(mapvar,2),'xticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                        caxset1=caxis;
                        caxset1(2)=max(max(mapvar))*0.9;
                        try
                        caxis(caxset1)
                        catch
                        end
                        title (['Laminar Amp ' ilabels{icik} ', ' num2str(wraw.avg.sweepno{trigcik}(1)) ' epochs, scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                        
                        ic1 = ic1+1;
                        subplot(imagecount2,4,ic1)
                        mapvar=mapvars_itc2_eeg{trigcik,icik};
                        
                        surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                        set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                        set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                        caxset1=caxis;
                        caxset1(2)=max(max(mapvar))*0.9;
                        caxis(caxset1)
                        title (['ITC ' ilabels{icik} ', ' num2str(wraw.avg.sweepno{trigcik}(1)) ' epochs, sig:' num2str(rayleigh_p(wraw.avg.sweepno{trigcik}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                        
                        ic1 = ic1+1;
                        subplot(imagecount2,4,ic1)
                        mapvar=mapvars_itc1_eeg{trigcik,icik};
                        surface(1:size(mapvar,2),1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                        set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                        set(gca,'xlim',[1 size(mapvar,2)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out','ydir','reverse')
                        set(gca,'xlim',[1 length(wraw.frq)],'xtick',1:ylabelres:size(mapvar,2),'xticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                        caxset1=caxis;
                        caxset1(2)=max(max(mapvar))*0.9;
                        caxis(caxset1)
                        title (['Laminar ITC ' ilabels{icik} ', ' num2str(wraw.avg.sweepno{trigcik}(1)) ' epochs, scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                        
                    end
                    
                    fname = filenamesout2{filecik}(1:end-4);
                    axes('Position',[0 0.98 1 0.2],'Visible','off');
                    text(0.5,0,[fname ', trigger locked averages, reprate = ' num2str(reprate(1)) ' Hz'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                    
                    print ('-djpeg', '-r1200', [directory4 'w02' num2str(trigcik) '_' filenamesout2{filecik}(1:end-4) '.jpg']);
                    
                end
                
            end
            
            close all
            
        end
        
        
        if ~isempty(wlfp.adrate)
            
            ylabels = num2str(wlfp.frq','%11.3g');
            time    = wlfp.time;
            
            imagecount_lfp = 3;
            if ~isempty(wlfp.avg.avg)
                imagecount_lfp2 = 3;
            else
                imagecount_lfp2 = 0;
            end
            
            curfig = figure;
            set(curfig,'position',[100   50   1200   800],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
            colormap jet
            ylabelres=10;
            fsize = 5;
            
            
            ic1 = 0;
            
            
            ic1=ic1+1;
            axpos=axposes1{2}(ic1,:);
            figureax1a=axes('Position',axpos);
            
            mapvar = squeeze(mean(wlfp.cnt_po(1:3,:,:),3));
            surface(1:size(mapvar,2),1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
            set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
            set(gca,'xlim',[1 size(mapvar,2)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out','ydir','reverse')
            set(gca,'xlim',[1 length(wlfp.frq)],'xtick',1:ylabelres:size(mapvar,2),'xticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
            
            cax1a=caxis;
            title([ 'Surface spectrogram,   cax: ' num2str(cax1a(1)) ' - ' num2str(cax1a(2))],'FontSize',8)
            
            ic1=ic1+1;
            axpos=axposes1{2}(ic1,:);
            figureax1a=axes('Position',axpos);
            
            mapvar =  squeeze(mean(wlfp.cnt_po(1:3,:,:),1));
            surface(time,1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
            set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
            set(gca,'xlim',[time(1) time(end)],'tickdir','out')
            %             if icik<imagecount
            %                 set(gca,'xticklabel','')
            %             end
            set(gca,'ylim',[1 length(wlfp.frq)],'ytick',1:ylabelres:size(mapvar,1),'yticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
            caxis(cax1a)
            cax3=caxis;
            title([ 'Surface time-frequency,   cax: ' num2str(cax3(1)) ' - ' num2str(cax3(2))],'FontSize',8)
            
            if ~isempty(wai.adrate)
                
                ylabels = num2str(wlfp.frq','%11.3g');
                time    = wai.time;
                
                ic1=ic1+1;
                axpos=axposes1{2}(ic1,:);
                figureax1a=axes('Position',axpos);
                
                mapvar = squeeze(mean(wai.cnt_po(2:4,:,:),3));
                surface(1:size(mapvar,2),1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
                set(gca,'xlim',[1 size(mapvar,2)],'ylim',[1 size(mapvar,1)],'ytick',1:1:size(mapvar,1),'tickdir','out','ydir','reverse')
                set(gca,'xlim',[1 length(wai.frq)],'xtick',1:ylabelres:size(mapvar,2),'xticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                
                cax1a=caxis;
                title([ 'Analog input spectrogram,   cax: ' num2str(cax1a(1)) ' - ' num2str(cax1a(2))],'FontSize',8)
                
                ic1=ic1+1;
                axpos=axposes1{2}(ic1,:);
                figureax1a=axes('Position',axpos);
                
                mapvar =  squeeze(mean(wai.cnt_po(2:4,:,:),1));
                surface(time,1:size(mapvar,1),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0],'color',[1 1 1])
                set(gca,'xlim',[time(1) time(end)],'ylim',[1 size(mapvar,1)],'tickdir','out')
                %                 if icik<imagecount
                %                     set(gca,'xticklabel','')
                %                 end
                set(gca,'ylim',[1 length(wai.frq)],'ytick',1:ylabelres:size(mapvar,1),'yticklabel',ylabels(1:ylabelres:end,:),'fontsize',fsize)
                caxis(cax1a)
                cax3=caxis;
                title([ 'Analog input time-frequency,   cax: ' num2str(cax3(1)) ' - ' num2str(cax3(2))],'FontSize',8)
                
            end
            
            
            fname = filenamesout2{filecik}(1:end-4);
            axes('Position',[0 0.98 1 0.2],'Visible','off');
            text(0.5,0,[fname ' - LFPAI'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
            
            print ('-djpeg', '-r1200', [directory4 'w01lfpai' num2str(trigcik) '_' filenamesout2{filecik}(1:end-4) '.jpg']);
            close all
            
            if imagecount_lfp2>0
                
                curfig = figure;
                set(curfig,'position',[100   50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                
                t3 = wlfp.avg.time;
                frq = wlfp.frq;
                ylabels     = num2str(frq','%11.3g');
                colormap jet
                ylabelres=10;
                fsize = 5;
                
                
                
                ic1 = 0;
                for icik = 1:imagecount_lfp2
                    
                    
                    ic1 = ic1+1;
                    subplot(imagecount_lfp2,2,ic1)
                    mapvar=squeeze(mean(wlfp.avg.po{trigcik}(1,icik,:,:),2));
                    surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                    set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                    set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                    caxset1=caxis;
                    caxset1(2)=max(max(mapvar))*0.9;
                    caxis(caxset1)
                    title (['Amp Surface-' num2str(icik) ', ' num2str(wlfp.avg.sweepno{trigcik}(1)) ' epochs, scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                    
                    ic1 = ic1+1;
                    subplot(imagecount_lfp2,2,ic1)
                    mapvar=squeeze(mean(wlfp.avg.itc{trigcik}(1,icik,:,:),2));
                    surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                    set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                    set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                    caxset1=caxis;
                    caxset1(2)=max(max(mapvar))*0.9;
                    if caxset1(1)>=caxset1(2)
                        caxset1(2)=caxset1(2)*1.1;
                    end
  
                    caxis(caxset1)
                    title (['ITC Surface-' num2str(icik) ', ' num2str(wlfp.avg.sweepno{trigcik}(1)) ' epochs, sig:' num2str(rayleigh_p(wlfp.avg.sweepno{trigcik}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                end
                fname = filenames{filecik}(1:end-4);
                axes('Position',[0 0.98 1 0.2],'Visible','off');
                text(0.5,0,[fname ' - Surf, trigger locked averages, reprate = ' num2str(reprate) ' Hz'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                
                print ('-djpeg', '-r1200', [directory4 'w02lfp' num2str(trigcik) '_' filenamesout2{filecik}(1:end-4) '.jpg']);
                
            end
            
            if ~isempty(wai.adrate) & imagecount_lfp2>0
                
                curfig = figure;
                set(curfig,'position',[100   50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto')
                
                t3 = wai.avg.time;
                frq = wai.frq;
                ylabels     = num2str(frq','%11.3g');
                colormap jet
                ylabelres=10;
                fsize = 5;
                
                
                
                ic1 = 0;
                for icik = 2:imagecount_lfp2+1
                    
                    
                    ic1 = ic1+1;
                    subplot(imagecount_lfp2,2,ic1)
                    mapvar=squeeze(mean(wai.avg.po{trigcik}(1,icik,:,:),2));
                    surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                    set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                    set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                    caxset1=caxis;
                    caxset1(2)=max(max(mapvar))*0.9;
                    caxis(caxset1)
                    title (['Amp Analog Input-' num2str(icik) ', ' num2str(wai.avg.sweepno{trigcik}(1)) ' epochs, scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                    
                    ic1 = ic1+1;
                    subplot(imagecount_lfp2,2,ic1)
                    mapvar=squeeze(mean(wai.avg.itc{trigcik}(1,icik,:,:),2));
                    surface(t3,1:size(frq,2),zeros(size(mapvar)),mapvar,'EdgeColor','none','FaceColor','interp');
                    set(gca,'xcolor',[.2 .2 .2],'ycolor',[.3 .3 .3],'color',[1 1 1],'fontsize',fsize)
                    set(gca,'xlim',[epoch_tframe(1) epoch_tframe(2)],'ylim',[1 size(frq,2)],'ytick',1:ylabelres:size(frq,2),'yticklabel',ylabels(1:ylabelres:end,:),'tickdir','out')
                    caxset1=caxis;
                    caxset1(2)=max(max(mapvar))*0.9;
                    try
                    caxis(caxset1)
                    catch
                    end
                    title (['ITC Analog Input-' num2str(icik) ', ' num2str(wai.avg.sweepno{trigcik}(1)) ' epochs, sig:' num2str(rayleigh_p(wai.avg.sweepno{trigcik}(1),0.05),'%11.2g') ', scale: ' num2str(caxset1(1),'%11.3g') '-' num2str(caxset1(2),'%11.3g')])
                end
                fname = filenames{filecik}(1:end-4);
                axes('Position',[0 0.98 1 0.2],'Visible','off');
                text(0.5,0,[fname ' - AI, trigger locked averages, reprate = ' num2str(reprate) ' Hz'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')
                
                print ('-djpeg', '-r1200', [directory4 'w02ai' num2str(trigcik) '_' filenamesout2{filecik}(1:end-4) '.jpg']);
                
            end
            
        end
        
        
        close all
        
    end
    
end

