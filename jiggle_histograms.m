figure
subplot(3,2,1)
histogram(DeltaMtxCoeffs((DeltaMtxCoeffs(:,4)<0.12),4),'BinWidth',0.02)
title('Jitter Tolerance MGB matrix')
hold on
subplot(3,2,3)
histogram(ThetaMtxCoeffs(ThetaMtxCoeffs(:,4)<0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,5)
histogram(AlphaMtxCoeffs(AlphaMtxCoeffs(:,4)<0.12,4),'BinWidth',0.02)
xlabel('Jitter Tolerance (Sec.)')
hold on
subplot(3,2,2)
histogram(DeltaMtxCoeffs((DeltaMtxCoeffs(:,4)==0.12),4),'BinWidth',0.02)
title('Tolerance > 0.12 S')
hold on
subplot(3,2,4)
histogram(ThetaMtxCoeffs(ThetaMtxCoeffs(:,4)==0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,6)
histogram(AlphaMtxCoeffs(AlphaMtxCoeffs(:,4)==0.12,4),'BinWidth',0.02)

figure
subplot(3,2,1)
histogram(DeltaOnCoeffs((DeltaOnCoeffs(:,4)<0.12),4),'BinWidth',0.02)
title('Jitter Tolerance across laminae (on BF)')
hold on
subplot(3,2,3)
histogram(ThetaOnCoeffs(ThetaOnCoeffs(:,4)<0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,5)
histogram(AlphaOnCoeffs(AlphaOnCoeffs(:,4)<0.12,4),'BinWidth',0.02)
xlabel('Jitter Tolerance (Sec.)')
hold on
subplot(3,2,2)
histogram(DeltaOnCoeffs((DeltaOnCoeffs(:,4)==0.12),4),'BinWidth',0.02)
title('Tolerance > 0.12 S')
hold on
subplot(3,2,4)
histogram(ThetaOnCoeffs(ThetaOnCoeffs(:,4)==0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,6)
histogram(AlphaOnCoeffs(AlphaOnCoeffs(:,4)==0.12,4),'BinWidth',0.02)



subplot(3,2,2)
histogram(DeltaOnCoeffs(DeltaOnCoeffs(:,4)<0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,4)
histogram(ThetaOnCoeffs(ThetaOnCoeffs(:,4)<0.12,4),'BinWidth',0.02)
hold on
subplot(3,2,6)
histogram(AlphaOnCoeffs(AlphaOnCoeffs(:,4)<0.12,4),'BinWidth',0.02)
legend('Delta','Theta','Alpha')
xlabel('Jitter Tolerance (Sec.)')
title('Jitter Tolerance across laminae (on BF)')