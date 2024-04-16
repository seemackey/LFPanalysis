%%%%%%%% important: bool.eegc = 1 only saves power, bool.eegc = 2 saves
%%%%%%%% phase this will be corrected


%%%%%%%% 2020 saccades03 hires-%% Noelle changed it from 5 3 3 to 4 3 3
close all; clear;tic
dir_nowave        = ['E:\spectrolaminar\AttnData\belt\cont\'];
dir_wave        = ['E:\spectrolaminar\AttnData\belt\wavelet\'];
filenames         = {};
wmethod           = [5 3 4]; %%% 9 7 3 high frq, 4 1 3 low, 5 1 3 all
frqrange          = [0.5 256];
epoch_tframe_orig = [-300 300];
bool.notch60hz    = 1;
bool.lfp          = 1;
bool.csd          = 1;
bool.mua          = 1;
bool.bipolar      = 0;
bool.hgamma       = 0;
bool.arej         = 0;
bool.eegc          = 0; 
bool.eegm          = 0; 
bool.eegb          = 0;
bool.eege          = 0;
bool.eegl          = 0; 
bool.eega          = 0; 
bool.eegh          = 0;
bool.downsample   = 0;
bool.ph           = 1;
bool.saccade_arej = 0;
bool.ttypes       = 1; %%%%%% if 0, all trigger types will be treated as the same,
%%%%%% if 1, they will be treated as individuals,
%%%%%% if 2, they will be renumbered based on distance from the deviant
%%%%%% if 3, epochs will be aligned to pattern onsets
bool.images = 1;
filestart = 1; % start with this file in your directory

[filenamesout2,filenames_nowave] =  snr_wavelet08cm(dir_nowave,dir_wave,filenames, wmethod, frqrange,  epoch_tframe_orig, bool,filestart);
beep on; beep;
toc