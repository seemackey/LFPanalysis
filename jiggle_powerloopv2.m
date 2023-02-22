%% loops through files and gets power pre stimulus and during stimulus

tic
clear
close all

%% set directory

mydir = 'H:\jitter\gt038_02\'; %gets directory
%mydir  = 'G:\jiggle\contproc\missing\gt044\wavelet\';
%mydir = 'F:\dyneyep\jitter\contproc\wavelet';
myfiles = dir(fullfile(mydir,'*@osw*')); %gets all files in struct
gran = [13]; % selecting gran channel, can be left blank for auto-select

%% loop through directory
for loopct = 1:length(myfiles)
    
    % file name stuff
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    
    % load wavelet file
    load(fullfilename);
   
    
    %% some inputs that determine which data we're going to use
 
    % found a block with two triggers
    if length(wraw.reprate)>1
        disp('OH NO! multiple stimulus streams')
    end
    
    % freq bands of interest ... hard coded for now, sorry, dear reader
    if wraw.reprate > 1 && wraw.reprate < 5 % delta
        
        frqbnd = 13;
        
    elseif wraw.reprate > 5.1 && wraw.reprate < 8 % theta
        
        frqbnd = 36;
        
    elseif wraw.reprate > 9 && wraw.reprate < 13 % alpha
        
        frqbnd = 45;
        
    end
    
         % pick time point of interest, median gets time zero (e.g. just before
        % stim onset
        time =  wraw.avg.time;
        mid = median(time);
        mididx = time(1,:)==mid;

        % find data before first trigger so we can calc baseline
        firsttrig = min(wraw.avg.trig{1,1}{1,1});
        pre = mean(wraw.cntc_po(:,frqbnd,1:round(firsttrig-0.01*firsttrig)),3);
        stim = wraw.avg.poc{1}(1,:,frqbnd,mididx)'; % power during stim
        
    if loopct == 1
        % diff between pre and stim
        powerdiff = stim-pre;

        % select channels based on diff
        supra = find(powerdiff==max(powerdiff(1:9)));
        if isempty(gran)
            gran = find(powerdiff==max(powerdiff(10:15)));
        end
        infra = find(powerdiff==max(powerdiff(17:21)));
    end
    
    stim_curr = wraw.avg.poc{1}(1,:,frqbnd,mididx)'; % power during stim
    % save supra, gran, and infra power during pre and stim
    output(loopct,:) = [pre(supra), pre(gran), pre(infra), stim_curr(supra), stim_curr(gran), stim_curr(infra)];

    clear wraw wlfp wai
end

% stick the channels we selected at the end so we can keep track of them
output(loopct+1,:) = [supra gran infra supra gran infra];

%% save important stuff
save(basefilename(1:12));
toc

%lets you know it's done
beep on; beep;

