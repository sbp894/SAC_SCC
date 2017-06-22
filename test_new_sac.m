clear;
close all;
clc;

addpath('mexsource\');

nReps=20;
DELAYbinwidth=50e-6;
rateMax=100;
nBoots=2e3;

[stim, fs]=audioread('dan_sent1.wav');
stimDur=length(stim)/fs;
Duration=.8*stimDur;

T(nBoots)=struct('old_sacc',nan, 'new_sacc', nan);

%%
CF_Hz=1e3;
ANmodel_Fs_Hz=100e3;
Cohc=1;
Cihc=1;
species=2;
SRtype_A=2;
noiseType=0;
implnt=0;

vIHC = model_IHC(stim.',CF_Hz,1,1/ANmodel_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
[~,meanrate_unad, ~] = model_Synapse(vIHC,CF_Hz,1,1/ANmodel_Fs_Hz,SRtype_A,noiseType ,implnt);
SpikeTrains=get_sptimes(meanrate_unad, ANmodel_Fs_Hz, nReps);
[temp,~,~,~] = ShufAutoCorr(SpikeTrains,DELAYbinwidth,Duration);
NSCCold=nan(nBoots, length(temp));
NSCCnew=nan(nBoots, length(temp));
%%
for montecarlo_loop=1:nBoots
    SpikeTrains=get_sptimes(meanrate_unad, ANmodel_Fs_Hz, nReps);
    
    tic;
    [NSCCold(montecarlo_loop, :),NSCCdelays_usec,AVGrates,TOTALspikes] = ShufAutoCorr(SpikeTrains,DELAYbinwidth,Duration);
    T(montecarlo_loop).old_sacc=toc;

    tic;
    [NSCCnew(montecarlo_loop, :),NSCCdelays_usec2,AVGrates2,TOTALspikes2] = SAChalf_m(SpikeTrains,DELAYbinwidth,Duration);
    T(montecarlo_loop).new_sacc=toc;
end

figure(1);
clf;
histogram([T.old_sacc],100);
hold on;
histogram([T.new_sacc],100);
legend('old', 'new');
title(sprintf('old mean time --> %1.4f, new mean time --> %1.4f', mean([T.old_sacc]), mean([T.new_sacc])));


figure(2);
plot(mean(NSCCold,1))
hold on;
plot(mean(NSCCnew,1), '--')
legend('old', 'new');