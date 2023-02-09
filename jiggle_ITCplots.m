%% plot ITC data

deltaon = [0.602408726
0.573605777
0.565511654
0.35069956
];

stdxax = [0,0.04,0.08,0.12];

deltaoff = [0.308057432
0.232273735
0.171754302
0.06662313

];

thetaon = [0.242877771
0.441017059
0.259340445
0.333335074

];

thetaoff = [0.238493178
0.127062021
0.203111992
0.022372152

];

alphaon = [0.795490559
0.366127121
0.176847849
0.18884624

];

alphaoff = [0.152629801
0.110657185
0.152591502
0.10288724

];
figure
subplot(1,3,1)
plot(stdxax,deltaon)
hold on
plot(stdxax,deltaoff)
legend('BF', 'non-BF')
ylabel('Intertrial Coherence')
xlabel('Std. Dev of Stim. Timing')
title('Delta ITC')

subplot(1,3,2)
plot(stdxax,thetaon)
hold on
plot(stdxax,thetaoff)
legend('BF', 'non-BF')
ylabel('Intertrial Coherence')
xlabel('Std. Dev of Stim. Timing')
title('Theta ITC')

subplot(1,3,3)
plot(stdxax,alphaon)
hold on
plot(stdxax,alphaoff)
legend('BF', 'non-BF')
ylabel('Intertrial Coherence')
xlabel('Std. Dev of Stim. Timing')
title('Alpha ITC')

