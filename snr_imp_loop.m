close all; clear all; 
tic; 

directory1 = ['/Volumes/Samsung03/pt019020/'];
  directory2 = ['/Volumes/16tbxfat/pt019020/'];
  filenames0   = {};
  chansets    = [1,23 ; 33,55];
%chansets    = [33,55];
  bool_preamp = 1;
   bool_avg = 1;

  [filenamesout] = snr_imp05(directory1,directory2,filenames0,chansets,bool_preamp, bool_avg);

  toc