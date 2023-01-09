function [eegcb, eegbb, trig01, tonefreqs, eegmb] = EphysExtractFxn(path)

%% extraction of lfp, csd, mua
% this function takes a file path and gets LFP and MUA data out of the file
% edits by Chase Mackey 2022


% load data
%path = ('C:\Users\cmackey\Documents\MATLAB\1-gt006000010@os_eye06_30');
load(path);
fname = path;


%% INPUTS (filtering, epoch timing etc.)%
timeframe_baseline              = [-200 0];
epoch_tframe            = [-250 250];  
newadrate =500;
filtere = [0.5 300];%LFP
filteru = [300 5000];%MUA
filtertype=1;
fsize = 6;
xlabelres=5;
%

%% epoching and filtering

% set timeframe
time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);

% extract continuous ephys data and triggers
[cnt_arej, cnte, cntm, cntc, cntu, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);
[trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig, params);

% bbn is trigger in the first position, led is trigger in second position
if isempty(trig.anatrig{1})
    trig0=trig.anatrig{2};
else
    trig0=trig.anatrig{1};
end


for trigredxct=1:length(trig0)
 trig01(trigredxct)    = round(trig0(trigredxct)./(trig.adrate/newadrate));
end


%filter the cnt in 1-100Hz to make nice CSD picts
%to bandpass filter in the 1.2-2.5Hz range the MUA (cntm)
n = 2;
Wn = [1 100]/(newadrate/2); %the first 2 numbers are the freq range you want to pass, which is =/- 20% of the reprate
[b,a] = butter(n,Wn);

% preallocate for filtered data
cntc_ff = zeros(size(cnte));
cntm_ff = zeros(size(cnte));
cnte_ff = zeros(size(cnte));
cntb_ff = zeros(size(cnte));

i=1;
for i=1:size(cntc,1) %this is for the electrode chs
 cntc_ff(i,:)=filtfilt(b,a,cntc(i,:)); %creates a new variable which contains filtered data
 cntm_ff(i,:)=filtfilt(b,a,cntm(i,:));
 cnte_ff(i,:)=filtfilt(b,a,cnte(i,:));
 cntb_ff(i,:)=filtfilt(b,a,cntb(i,:));
end

%epoching
x1 = round(epoch_tframe(1)*(newadrate/1000));
x2 = round(epoch_tframe(2)*(newadrate/1000));

i=1;
for i=1:length(trig01)
    
 eegm(:,i,:)=cntm_ff(:,trig01(i)+x1:trig01(i)+x2);%epoching the mua
 eegc(:,i,:)=cntc_ff(:,trig01(i)+x1:trig01(i)+x2);%epoching the csd
 eege(:,i,:)=cnte_ff(:,trig01(i)+x1:trig01(i)+x2);% same for lfp
 eegb(:,i,:)=cntb_ff(:,trig01(i)+x1:trig01(i)+x2);% same for bipolar lfp
 
end


%% baseline correct
chct=1; % channels
trct=1; % trials

for chct=1:size(eegm,1)
 for trct=1:size(eegm,2)
     eegmb(chct,trct,:)   = squeeze(eegm(chct,trct,:))-squeeze(mean(eegm(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegcb(chct,trct,:)   = squeeze(eegc(chct,trct,:))-squeeze(mean(eegc(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegeb(chct,trct,:)   = squeeze(eege(chct,trct,:))-squeeze(mean(eege(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
     eegbb(chct,trct,:)   = squeeze(eegb(chct,trct,:))-squeeze(mean(eegb(chct,trct,max(find(time<=timeframe_baseline(1))):max(find(time<=timeframe_baseline(2)))),3));
 end
end


%triggers if there are multiple
if isempty(trig.ttype(1))
    tonefreqs=trig.ttype(2);
else
    tonefreqs=trig.ttype(1);
end

end