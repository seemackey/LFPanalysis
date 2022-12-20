function [eegcb, eegbb, trig01, OutlierIdxUnique] = rejectartifacts(eegcb,eegbb,trig01)
    %% Artifact reject
    % takes baseline corrected and epoched lfp, csd, mua and returns clean data
    %preallocate
    %mmax = zeros(size(eegm,2),1);
    cmax = zeros(size(eegcb,2),1);
    bmax = zeros(size(eegbb,2),1);

    for i1=1:size(eegcb,2)

     %mmax(i1)=max(mean(abs(squeeze(eegmb(:,i1,:))),1));
     cmax(i1)=max(mean(abs(squeeze(eegcb(:,i1,:))),1));
     bmax(i1)=max(mean(abs(squeeze(eegbb(:,i1,:))),1));

    end

    %[b,OutlieridxM,outliers] = deleteoutliers(mmax,0.1);%idXM are the trigs with artifacts
    [b,OutlieridxC,outliers] = deleteoutliers(cmax,0.1);%
    [b,OutlieridxB,outliers] = deleteoutliers(bmax,0.1);%
    OutlierIdx=[OutlieridxC; OutlieridxB];
    OutlierIdxUnique = unique(OutlierIdx);

    %% outputs of fxn. baseline corrected, outliers removed etc. 
    trig01(:,OutlierIdxUnique) = []; %trigger
    %eegmb(:,OutlierIdxUnique,:) = []; %  MUA 
    eegcb(:,OutlierIdxUnique,:) = []; %  CSD 
    eegbb(:,OutlierIdxUnique,:) = []; %  bip LFP 

end