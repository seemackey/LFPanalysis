 function [filenamesout] = snr_imp05(directory1,directory2,filenames0,chansets,bool_preamp, bool_avg)

% close all; clear all; tic; directory1 = ['G:\Original_Data\Robin\exps\rb017018\'];
%   directory2 = ['E:\Imported_data\Robin\contproc\rb017018\'];
%   filenames0   = {};
%   chansets    = [1,23 ; 33,55];
%   bool_preamp = 1;
%    bool_avg = 1;
%   [filenamesout] = snr_imp05(directory1,directory2,filenames0,chansets,bool_preamp, bool_avg);
%   toc


% clear all; close all; tic
% directory1 = 'G:\Original_Data\Gertie\exps\x\';
% dirnames=[];
% f=dir( [ directory1]);
% for i=1:1:size(f,1)
%     ff=f(i).name;
%     dirnames{i}=ff;
% end
% dirnames=dirnames(3:end);
% dirnames = dirnames';
% for dircik = 1:length(dirnames);
%     tic; directory01  = [directory1 dirnames{dircik} '\'];
%     directory02  = ['E:\Imported_data\Gertie\contproc\' dirnames{dircik} '\'];
%     filenames0   = {};
%     chansets    = [1,23 ; 33,55];
%     bool_preamp = 1;
%     bool_avg = 1;
%     [filenamesout] = snr_imp05(directory01,directory02,filenames0,chansets,bool_preamp, bool_avg);
%     toc
% end


if isempty(filenames0)
    filenames=[];
    f=dir( [ directory1 '*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
else
    filenames=[];
    f=dir( [ directory1 filenames0{1} '*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        filenames{i}=ff;
    end
end

i1=1;
filenames2=filenames;
while i1<length(filenames2)
    if filenames2{i1}(1:11)==filenames2{i1+1}(1:11)
        filenames2(i1+1)=[];
    else
        i1=i1+1;
    end
end

filecount=0;

for maincik=1:size(filenames2,2)
    if maincik==1
        [~,~,~] = mkdir(directory2);
    end
    
    
    
    cnt0        = [];
    craw.cnt    = [];
    cai.cnt     = [];
    clfp.cnt    = [];
    craw.adrate = [];
    cai.adrate  = [];
    clfp.adrate = [];
    
    fns=[];
    fns_event=[];
    f=dir( [ directory1 filenames2{maincik}(1:11) '*.mat']);
    for i=1:1:size(f,1)
        ff=f(i).name;
        fns{i}=ff;
    end
    
    for i1=1:8
        trig.anatrig{i1}       = [];
        trig.triglength{i1}    = [];
        trig.ttype{i1}         = [];
    end
    trig.digtrig            = [];
    otherdata               = [];
    params.oldnew           = 2;
    trig.trigchan           = [];
    
    
    %%%% concatenate triggers and recordings across file segments
    for filecik = 1:length(fns)
        
        load([directory1 fns{filecik}])
        all_loaded_variables = whos('-file',[directory1 fns{filecik}]);
       
        
        %%% trigger
        if filecik == 1
            try
                fbegin = CRAW_001_TimeBegin;
            catch
                fbegin = CAI_001_TimeBegin;
            end
        end
        
        for i1=1:4
            te1=['CTTL_' num2str(i1,'%03.0f')];
            if exist([te1 '_TimeBegin'],'var')==1
                trig.adrate=eval([te1 '_KHz'])*1000;
                tr1 = eval([te1 '_Up']);
                tr2 = eval([te1 '_Down']);

                if length(tr1)>length(tr2)
                    tr1=tr1(1:length(tr2));
                elseif length(tr2)>length(tr1)
                    tr2=tr2(1:length(tr1));
                end
                trig.anatrig{i1}=[trig.anatrig{i1} (eval([te1 '_TimeBegin'])-fbegin)*trig.adrate+tr1];
                trig.triglength{i1}=[trig.triglength{i1} (tr2-tr1)];
                
            end
        end
        
        for i1=1:4
            
            te1=['CStimMarker_' num2str(i1,'%03.0f')];
            if exist([te1],'var')==1
                trig.adrate=eval([te1 '_KHz'])*1000;
                tr1 = eval(te1);
                trig.anatrig{i1+4}=[trig.anatrig{i1+4} (tr1(1,:)/trig.adrate-fbegin).*trig.adrate];
                trig.triglength{i1+4}=[trig.triglength{i1+4} tr1(2,:)];
            end
        end
       
        
        
        %%% cnt, if chansets is not empty
        if ~isempty(chansets)
            
            b=[];
            for i1=1:size(chansets,1)
                chansall(i1,:) = [length(b)+1 length(b)+length(chansets(i1,1):chansets(i1,2))];
                b=[b chansets(i1,1):chansets(i1,2)];
            end
            
            craw.adrate = eval(['CRAW_' num2str(b(1),'%03.0f') '_KHz'])*1000;
            cnt=zeros(length(b),length(CRAW_001));
            for i1=1:length(b)
                cnt(i1,:)=[double(eval(['CRAW_' num2str(b(i1),'%03.0f')]))];
            end
            cnt0=[cnt0,cnt];
            clear cnt
        end
        
        
        
        %%% analog input
        try
            b = 1:5;
            cai.adrate = eval(['CAI_' num2str(b(1),'%03.0f') '_KHz'])*1000;
            cnt=zeros(length(b),length(CAI_001));
            for i1=1:length(b)
                cnt(i1,:)=[double(eval(['CAI_' num2str(b(i1),'%03.0f')]))];
            end
        catch
            
        end
        if exist ('cnt','var')
            cai.cnt=[cai.cnt,cnt];
            clear cnt
        end
        
        
        %%% EEG
        try
            b = 65:69;
            clfp.adrate = eval(['CLFP_' num2str(b(1),'%03.0f') '_KHz'])*1000;
            cnt=zeros(length(b),length(CLFP_065));
            for i1=1:length(b)
                cnt(i1,:)=[double(eval(['CLFP_' num2str(b(i1),'%03.0f')]))];
            end
            
        end
        
        try
            b = 1:5;
            clfp.adrate = eval(['CLFP_' num2str(b(1),'%03.0f') '_KHz'])*1000;
            cnt=zeros(length(b),length(CLFP_001));
            for i1=1:length(b)
                cnt(i1,:)=[double(eval(['CLFP_' num2str(b(i1),'%03.0f')]))];
            end
            
        end
        if exist ('cnt','var')
            clfp.cnt=[clfp.cnt,cnt];
            clear cnt
        end
        
        clear (all_loaded_variables.name)
    end
    
    
    for i1=1:length(trig.anatrig)
        if ~isempty(trig.anatrig{i1})
            trig.trigchan=[trig.trigchan i1];
        end
    end
    
    try
        
        x=importdata([ directory1 filenames2{maincik}(1:end-9) '@v.dat']);
        y=importdata([ directory1 filenames2{maincik}(1:end-9) '@y.dat']);
        
        fns_event{1}=[filenames2{maincik}(1:end-9) '@v.dat'];
        fns_event{2}=[filenames2{maincik}(1:end-9) '@y.dat'];
       
        
        if size(x,2)>length(trig.trigchan) && size(x,2)<15
            ciklength = length(trig.trigchan);
            params.filedat=y.data;
            params.filedat_label=y.textdata;
        elseif size(x,2)<=length(trig.trigchan)
            ciklength = length(trig.trigchan);
            params.filedat=y.data;
            params.filedat_label=y.textdata;
        else
            ciklength = 0;
            if size(x,2)==66
                trig.ttype = x(:,1:16);
                trig.digtrig = x(:,17:32)./100000*trig.adrate;
                params.loud1 = x(:,33:48);
                params.loud2 = x(:,49:64);
                params.buttpress = x(:,65:66);
            elseif size(x,2)==63
                trig.ttype = [x(:,1:16)];
                trig.digtrig = [x(:,17:31),zeros(size(x,1),1)]./100000*trig.adrate;
                params.loud1 = x(:,32:46);
                params.loud2 = x(:,47:61);
                params.buttpress = x(:,62:63);
            end
            params.filedat=y.data;
            params.filedat_label=y.textdata;
            otherdata.eventfile=2;
        end
        
        for i=1:ciklength
            %             if i<=ciklength
            if  size(1:length(trig.anatrig{trig.trigchan(1)}),2)<size(x,1)
                trig.ttype{trig.trigchan(i)}=x(1:length(trig.anatrig{trig.trigchan(i)}),i);
            else
                trig.ttype{trig.trigchan(i)}=x(:,i);
                trig.ttype{trig.trigchan(i)}(end+1:length(trig.anatrig{trig.trigchan(i)}))=0;
            end
            %             else
            %                 trig.ttype{i}=[];
            %                 trig.triglength{i}=[];
            %             end
        end
        otherdata.eventfile=1;
        
    catch
        
        for i=1:size(trig.trigchan,2)
            trig.ttype{trig.trigchan(i)}=ones(length(trig.anatrig{trig.trigchan(i)}),1);
        end
        otherdata.eventfile=0;
        
    end

    
    otherdata.chansets = chansets;
    
    otherdata.bool_preamp = bool_preamp;
    
    if ~isempty(chansets)
        
        if bool_preamp == 1
            cnt0 = cnt0./10;
        end
        
        for findex=1:size(chansets,1)
            filecount=filecount+1;
            chans=chansall(findex,:);
            if size(chansets,1)>1
                findex1=findex;
            else
                if chansets(1)>24 && chansets(1)<49
                    findex1=2;
                elseif chansets(1)>48
                    findex1=3;
                else
                    findex1=1;
                end
            end
            craw.cnt    = cnt0(chans(1):chans(2),:,:);
            
            filenamesout{filecount} = [num2str(findex1) '-' filenames2{maincik}(1:end-9) '@os.mat'];
            
            craw.fname = [num2str(findex1) '-' filenames2{maincik}(1:end-9)];
            if bool_avg == 1
                filtere = [0.5 300];
                filteru = [300 5000];
                epoch_tframe = [-25 250];
                image_tframe = [-25 250];
                bool_baseline = 1;
                newadrate = craw.adrate;
                newadrate = 20000;
                filtertype = 1;
                [craw,clfp,cai] = module_avg01(craw, clfp, cai, trig, params, newadrate, filtere, filteru, filtertype, epoch_tframe,  bool_baseline);
                timeframe = [-5 50];
                [craw] = module_avgfig01(craw,timeframe);
                
                print('-djpeg','-r300',[directory2 filenamesout{filecount}(1:end-4) '.jpg'])
                close all
            end
            
            
            otherdata.fns=fns;
            otherdata.fns_event=fns_event;
            
            
            
            save([directory2 filenamesout{filecount}], 'craw','cai','clfp','trig','otherdata','params','-mat','-v7.3')
        end
        
        clear cnt0
        
    else
        filecount=filecount+1;
        chansall = 100;
        if chansall(1)>24 && chansall(1)<49
            findex1=2;
        elseif chansall(1)>48
            findex1=0;
        else
            findex1=1;
        end
        
        craw.fname = [num2str(findex1) '-' filenames2{maincik}(1:end-9)];
        if bool_avg == 1;
            filtere = [0.5 300];
            filteru = [300 5000];
            epoch_tframe = [-25 250];
            image_tframe = [-5 50];
            bool_baseline = 1;
            newadrate = craw.adrate;
            filtertype = 1;
            [craw,clfp,cai] = module_avg01(craw, clfp, cai, trig, params, newadrate, filtere, filteru, filtertype, epoch_tframe,  bool_baseline);
        end
        
        
        filenamesout{filecount} = [num2str(findex1) '-' filenames2{maincik}(1:end-9) '@os.mat'];
        
        otherdata.fns=fns;
        otherdata.fns_event=fns_event;
        
        save([directory2 filenamesout{filecount}], 'craw','cai','clfp','trig','otherdata','params','-mat','-v7.3')
        
    end
    
    
    clear cnt0 craw cai clfp trig
    close all
end

