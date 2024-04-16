% filter for unit before spike sorting
clear
clc;
path = ['C:\Users\cmackey\Documents\AttnTuning\ab016\2-ab015016018@sortedv2.mat'];

load(path);

newadrate = 9000;
filtere = [0.5 300];%LFP
filteru = [300 3000];%MUA
filtertype=1;
trigch=1;

[cnt_arej, ~, ~, ~, cntu, ~] = module_cnt05cm(craw, newadrate, filtere, filteru, filtertype);

% Modify the path to replace "@os" with "@su"
newPath = strrep(path, '@os', '@su');

% Save the cntu variable to the new path
save(newPath, 'cntu');

% Modify the path to replace "@su" with "@trigs"
trigPath = strrep(newPath, '@su', '@trigs');

% Save the cntu variable to the new path
save(trigPath, 'trig');

%[trigtype,trigtimes] = EphysExtractTrigs(trig,trigch,newadrate);
% function [trig1s,ttype1s, triglength1s, findex2s] = module_trig01(trig,params);
