%% calculating trigger onset asynchrony

%% path here
%path = '';
%load(path);

%% 
presenttrigs = trig.digtrig>0;
ttimes = trig.digtrig(presenttrigs);
ttimes = ttimes*2; % convert to ms, AD rate was 500

for trigct = 1:1:length(ttimes)
    
    soa(trigct) = ttimes(trigct)-((1000/wraw.reprate)*(trigct));
    
end

histogram(soa,'BinWidth',30)
meansoa = mean(soa)
sdsoa = std(soa)