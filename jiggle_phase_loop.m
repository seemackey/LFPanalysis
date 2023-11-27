% This loops through a directory of wavelet files and extracts ITC
tic
clear
close all

%% set directory
% tic
mydir = '/Volumes/14TB_USB_01/jitter/rb037_02/'; %gets directory
%mydir  = 'G:\jiggle\contproc\missing\gt044\wavelet\';
%mydir = 'F:\dyneyep\jitter\contproc\wavelet';
myfiles = dir(fullfile(mydir,'*.mat')); %gets all files in struct
chans_sel = 5;
area = 2; % 1 for mgb, 2 for a1

%% loop through directory

for loopct = 1:length(myfiles)
    
    % file name stuff
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    
    % load wavelet file
    load(fullfilename);
   
    
    %% some inputs that determine which ITC values we're going to use
    % time re: stim onset
   
    time =  wraw.avg.time;
    
    
    rate = wraw.reprate(1);
    
    % found a couple blocks with two triggers
    if length(wraw.reprate)>1
        disp('OH NO! multiple triggers')
    end
    
    % freq bands of interest ... hard coded for now, sorry, dear reader
    % CHANGED FOR MGB so put an if statement here to deal wit hit
    if max(wraw.frqrange)>50 % if the max freq was high, used these
        if rate > 1 && rate < 5 % delta

            frqbnd = 13; % 13 or 21

        elseif rate > 5.1 && rate < 8 % theta

            frqbnd = 36; %36 or 43

        elseif rate > 9 && rate < 13 % alpha

            frqbnd = 45; %45 or 53

        end
    end
    
    if max(wraw.frqrange)<51 % if the max freq was low, used these
        if rate > 1 && rate < 5 % delta

        frqbnd = 21; % 13 or 21

        elseif rate > 5.1 && rate < 8 % theta

        frqbnd = 43; %36 or 43

        elseif rate > 9 && rate < 13 % alpha

        frqbnd = 53; %45 or 53

        end
    end
    
    
    %% epoch continuous signal in wraw file
    
    [eegbb, trig01, trigtype, eegmb,eegphase] = EphysEpochFxnWithPhase(wraw,trig,1,30,area);
    
    
    %% number of channels to loop through
    
    if area ==1
        numchans = length(wraw.avg.itcb{1}(1,:,1,1));
    elseif area ==2
        numchans = length(wraw.avg.itcc{1}(1,:,1,1));
    end
    %numchans = length(wraw.avg.itcm{1}(1,:,1,1));
    
    
    % pick time point of interest, median gets time zero (e.g. just before
    % stim onset
   mid = median(time);
   mididx = time(1,:)==mid;
    
    for chanct = 1:1:numchans
        
        phase.onetp{loopct,:,chanct} = eegphase(chanct,:,17); % itc at tp of interest
        
        phase.frqchosen = wraw.frq(frqbnd);
        
        mua{loopct,:,chanct} = eegmb(chanct,:,15:30);
    end
    
    [meansoa,sdsoa] = jiggle_soacalc(wraw,trig);
    phase.jitter(loopct) = sdsoa';
    
    %itc.selectch{loopct,1} = mean(squeeze(wraw.avg.itcc{1}(1,chans_sel,frqbnd,mididx))); % itc at tp of interest, avg across channels
    %itc.std{loopct,1} = std(squeeze(wraw.avg.itcc{1}(1,chans_sel,frqbnd,mididx))); % itc at tp of interest, std dv across channels
    
    clear wraw wlfp wai
end
toc

%% save important stuff

save(basefilename(1:12));
beep on; beep;