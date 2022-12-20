 function [trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig,params);

streamon_digi = find(sum(trig.digtrig,1)~=0);
try
streamon_ana  = trig.trigchan;
catch
    streamon_ana  = [];
end

a   = 0;
del = [];
for i1=1:length(streamon_ana)
    if length(trig.anatrig{streamon_ana(i1)})<20
        a=a+1;
        del(a)=i1;
    end
end


streamon_ana(del)=[];

% a   = 0;
% del = [];
% for i1=1:length(streamon_digi)
%     if length(trig.digtrig{i1})<5
%         a=a+1;
%         del(a)=i1;
%     end
% end
% 
% streamon_digi(del)=[];

if ~isempty(streamon_digi)  & isempty(trig.anatrig{3})
    streamon = streamon_digi;
    trigswitch = 1;
else
    streamon = streamon_ana;
    trigswitch = 0;
end

findex2s        = [];
trig1s          = [];
ttype1s         = [];
triglength1s    = [];


for trigcik = 1:length(streamon)
    
    if trigswitch == 0
        trig1           = trig.anatrig{streamon(trigcik)};

        try
            ttype1          = trig.ttype{streamon(trigcik)}(1:length(trig1));
        catch
            
                ttype1          = ones(1,length(trig1));


        end
        
        try
            triglength1     = trig.triglength{streamon(trigcik)};
            
        catch
               triglength1     = ones(1,length(trig1));
            end
      
        
        if length(triglength1)<length(trig1)
            triglength1=[triglength1,0];
        end
        
        
    elseif trigswitch == 1
        
        trig1=trig.digtrig(:,streamon(trigcik));
        ttype1=trig.ttype(:,streamon(trigcik));
        triglength1 = ones(size(trig1));
        if trig1(2)~=0
            x=trig1(2:end);
            x(x==0)=[];
            trig1=[trig1(1);x];
            
            x=ttype1(2:end);
            x(x==0)=[];
            ttype1=[ttype1(1);x];
            
        else
            trig1(trig1==0)=[];
            ttype1(ttype1==0)=[];
        end
        
        if ~isempty(trig.anatrig{1}) && ~isempty(trig.anatrig{2})
            x=[trig.anatrig{1}(1) trig.anatrig{2}(1)];
        elseif ~isempty(trig.anatrig{1})
            x=[trig.anatrig{1}(1)];
        elseif ~isempty(trig.anatrig{2})
            x=[trig.anatrig{2}(1)];
        end
        if streamon(trigcik)<16
            trig1=trig1+min(x);
        else
            trig1=trig1+min(x);
        end
        
        
    end
    
    if size(ttype1,2)>size(ttype1,1)
        ttype1=ttype1';
    end
    
    if length(ttype1)>length(trig1)
        ttype1=ttype1(1:length(trig1));
    end

    
    if streamon(trigcik)<16 && trigswitch == 1
        findex2=['_' params.filedat_label{1} '_' num2str(round(params.filedat(6,streamon(trigcik)))) '_' num2str(round(params.filedat(1,streamon(trigcik))*1000))];
    elseif streamon(trigcik)==16 && trigswitch == 1
        findex2=['_' params.filedat_label{1} '_' 'led' '_' num2str(round(params.filedat(1,streamon(trigcik))*1000))];
    else
        findex2=num2str(streamon(trigcik));
    end
    
    findex2s{trigcik}       = findex2;
    trig1s{trigcik}         = trig1;
    ttype1s{trigcik}        = ttype1;
    triglength1s{trigcik}   = triglength1;
    
end

%%%% create extra streams for patterns
try
    if strcmp(params.filedat_label{1} , 'SpeedTono')

        findex2s{1}=['_SpTono_all'];
        findex2s{2}=['_SpTono_pa1'];
        findex2s{3}=['_SpTono_pa2'];
        
        z=find(ttype1>1000);
        
      
        
        trigcik = 2;
        trig1s{trigcik}   = trig1(z(1:11:end));
        ttype1s{trigcik}  = ttype1(z(1:11:end));
        triglength1s{trigcik} = triglength1(z(1:11:end));
        
        trigcik = 3;
        trig1s{trigcik}    = trig1(z(1:11*5:end));
        ttype1s{trigcik}  = ttype1(z(1:11*5:end));
        triglength1s{trigcik} = triglength1(z(1:11*5:end));
        
    end
catch
end

%%%% create extra strreams for pre and pos-led stimuli
try
    if length(findex2s{1})>5
        if strcmp(findex2s{1}(1:6), '_NewAV')
            
            for i1=1:length(trig1s)
                x1(i1)=round(median(diff(trig1s{i1})));
            end
            
            paired_stream=find(x1==x1(end));
            
            if length(paired_stream)>1
                
                trig1=trig1s{paired_stream(1)};
                ttype1=ttype1s{paired_stream(1)};
                triglength1=triglength1s{paired_stream(1)};

                trig2=trig1s{paired_stream(2)};
                
                led_start=max(find(trig1==trig2(1)));
                led_end=max(find(trig1==trig2(end)));
                
                trigcik = length(trig1s)+1;
                trig1s{trigcik}    = trig1(1:led_start-1);
                ttype1s{trigcik}  = ttype1(1:led_start-1);
                triglength1s{trigcik} = triglength1(1:led_start-1);
                findex2s{trigcik}=[findex2s{paired_stream(1)} '-1'];
                
                trigcik = length(trig1s)+1;
                trig1s{trigcik}    = trig1(led_end+1:end);
                ttype1s{trigcik}  = ttype1(led_end+1:end);
                triglength1s{trigcik} = triglength1(led_end+1:end);
                findex2s{trigcik}=[findex2s{paired_stream(1)} '-2'];
                
            end
        end
    end
catch
end