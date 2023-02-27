%% looking at pupil response


%% load data
path = '/Users/chase/Desktop/NKI/data/1-bu037038018@os_eye06_20.mat';
[~, ~, trigtimes, trigtypes, ~,pupepdb,pupepd] = EphysAndPupilExtractFxn(path,1);



%% first check time-course of pupil resp to std vs deviant stimuli
pupil_avgtc = mean(pupepdb(1,:,:),1)'; % average MUA across chs and whole epoch
uniquetrigs = unique(trigtype);
    
    
    % sort by trigger type (e.g. tone frequency or led type
    for trig_ct = 1:1:length(uniquetrigs)
        
        currtrig = trigtype(:,1)==uniquetrigs(trig_ct,1);
        pupavgtcsort(trig_ct,1) = mean(pupil_avgtc(currtrig==1,1));
        pupstdtcsort(trig_ct,1) = std(pupil_avgtc(currtrig==1,1));
        
    end
    
    
%% now check "tuning" of pupil
pupil_avg = mean(mean(pupepdb(1,:,tpmua1:tpmua2),1),3)'; % average MUA across chs and whole epoch
uniquetrigs = unique(trigtype);
    
    
    % sort by trigger type (e.g. tone frequency or led type
    for trig_ct = 1:1:length(uniquetrigs)
        
        currtrig = trigtype(:,1)==uniquetrigs(trig_ct,1);
        MUAavgsort(trig_ct,1) = mean(MUAavg(currtrig==1,1));
        MUAstdsort(trig_ct,1) = std(MUAavg(currtrig==1,1));
        
    end