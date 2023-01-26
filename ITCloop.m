% This loops through the ephys extract fxn and gets MUA out, then checks 
% if the response is significantly greater than baseline

clear
close all

%% set directory
tic
mydir = 'C:\Users\cmackey\Documents\jiggle\wavelet'; %gets directory
myfiles = dir(fullfile(mydir,'*.mat')); %gets all files in struct


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
    
    % freq bands of interest ... hard coded for now, sorry, dear reader
    if wraw.reprate > 1 && wraw.reprate < 5 % delta
        
        frqbnd = 13;
        
    elseif wraw.reprate > 5.1 && wraw.reprate < 8 % theta
        
        frqbnd = 36;
        
    elseif wraw.reprate > 9 && wraw.reprate < 13 % alpha
        
        frqbnd = 45;
        
    end
    
    % number of channels to loop through
    numchans = length(wraw.avg.itcc{1}(1,:,1,1));
    
    % pick time point of interest, median gets time zero (e.g. just before
    % stim onset
   mid = median(time);
   mididx = time(1,:)==mid;
    
    for chanct = 1:1:numchans
        itc.all{loopct,chanct} = squeeze(wraw.avg.itcc{1}(1,chanct,frqbnd,:)); % itc across time
        itc.onetp{loopct,chanct} = squeeze(wraw.avg.itcc{1}(1,chanct,frqbnd,mididx)); % itc at tp of interest
        
    end
    itc.mean{loopct,1} = mean(squeeze(wraw.avg.itcc{1}(1,:,frqbnd,mididx))); % itc at tp of interest, avg across channels
    itc.std{loopct,1} = std(squeeze(wraw.avg.itcc{1}(1,:,frqbnd,mididx))); % itc at tp of interest, std dv across channels
end
toc

%% save important stuff
clear wraw wlfp wai

save('jigglesite1');