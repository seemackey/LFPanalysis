%% loop through itc values and fit with exp function
close all
makecurve=1;
innerct=1;
for fittingct = 1:4:length(alphaon)

     y=alphaon(fittingct:fittingct+3);
%      figure
    [cfs,jit_tol] = taufxn_v3(sd,y,0.01,0.001,makecurve);
    AlphaOnCoeffs(innerct,1:3) = cfs;
    AlphaOnCoeffs(innerct,4) = jit_tol(1);
    AlphaOnCoeffs(innerct,5) = jit_tol(2);
    innerct=innerct+1;
end

figure
subplot(1,2,1)
histogram(AlphaOnCoeffs(:,4),'BinWidth',0.01)
xlabel('Jitter Tolerance (Sec.)')
title('Jitter Tolerance across laminae (on BF)')
subplot(1,2,2)
histogram(AlphaOnCoeffs(:,5))
title('Channels without entrainment')


%% plotting


figure
subplot(1,2,1)
histogram(DeltaOffCoeffs(:,4),'BinWidth',0.01)
hold on
histogram(ThetaOffCoeffs(:,4),'BinWidth',0.01)
hold on
histogram(AlphaOffCoeffs(:,4),'BinWidth',0.01)
legend('Delta','Theta','Alpha')
xlabel('Jitter Tolerance (Sec.)')
title('Jitter Tolerance across laminae (off BF)')
subplot(1,2,2)
histogram(DeltaOnCoeffs(:,4),'BinWidth',0.01)
hold on
histogram(ThetaOnCoeffs(:,4),'BinWidth',0.01)
hold on
histogram(AlphaOnCoeffs(:,4),'BinWidth',0.01)
legend('Delta','Theta','Alpha')
xlabel('Jitter Tolerance (Sec.)')
title('Jitter Tolerance across laminae (on BF)')

figure
subplot(1,2,1)
cdfplot(DeltaOffCoeffs(:,4))
hold on
cdfplot(ThetaOffCoeffs(:,4))
hold on
cdfplot(AlphaOffCoeffs(:,4))
legend('Delta','Theta','Alpha')
xlabel('Jitter Tolerance (Sec.)')
ylabel('Cumulative Proportion')
title('Jitter Tolerance across laminae (off BF)')
grid off
subplot(1,2,2)
cdfplot(DeltaOnCoeffs(:,4))
hold on
cdfplot(ThetaOnCoeffs(:,4))
hold on
cdfplot(AlphaOnCoeffs(:,4))
legend('Delta','Theta','Alpha')
xlabel('Jitter Tolerance (Sec.)')
ylabel('Cumulative Proportion')
title('Jitter Tolerance across laminae (on BF)')
grid off
