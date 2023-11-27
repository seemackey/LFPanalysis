%% calculating trigger onset asynchrony

%% path here
%path = '';
%load(path);

%% 
function [meansoa,sdsoa] = jiggle_soacalc(wraw,trig)
    presenttrigs = trig.digtrig>0; 
    ttimes = trig.digtrig(presenttrigs);
    ttimes = ttimes*(1000/wraw.adrate);
    soadiff=diff(ttimes);
    soa=soadiff-(1000/wraw.reprate);


    %histogram(soa,'BinWidth',10);
    meansoa = mean(soa);
    sdsoa = std(soa);
end