clear;
close all;
clc;

addpath('mexsource\');

nReps=20;
DELAYbinwidth=50e-6;
nBoots=.3e3;

[stim, fs]=audioread('dan_sent1.wav');
stimDur=length(stim)/fs;
SCCdur=.8*stimDur;

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

vIHC_p= model_IHC(stim.',CF_Hz,1,1/ANmodel_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
[~,meanrate_unad_pos, ~] = model_Synapse(vIHC_p,CF_Hz,1,1/ANmodel_Fs_Hz,SRtype_A,noiseType ,implnt);

vIHC_n= model_IHC(-stim.',CF_Hz,1,1/ANmodel_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
[~,meanrate_unad_neg, ~] = model_Synapse(vIHC_p,CF_Hz,1,1/ANmodel_Fs_Hz,SRtype_A,noiseType ,implnt);


SpikeTrains_p=get_sptimes(meanrate_unad_pos, ANmodel_Fs_Hz, nReps);
SpikeTrains_n=get_sptimes(meanrate_unad_neg, ANmodel_Fs_Hz, nReps);
[temp,~,~,~] = ShufCrossCorr({SpikeTrains_p, SpikeTrains_n},DELAYbinwidth,SCCdur);
NSCCold=nan(nBoots, length(temp));
NSCCnew=nan(nBoots, length(temp));
%%
for montecarlo_loop=1:nBoots
    SpikeTrains_p=get_sptimes(meanrate_unad_pos, ANmodel_Fs_Hz, nReps);
    SpikeTrains_n=get_sptimes(meanrate_unad_neg, ANmodel_Fs_Hz, nReps);
    
    tic;
    [NSCCold(montecarlo_loop, :),NSCCdelays_usec,AVGrates,TOTALspikes] = ShufCrossCorr({SpikeTrains_p, SpikeTrains_n},DELAYbinwidth,SCCdur);
    T(montecarlo_loop).old_sacc=toc;
    
    tic;
    [NSCCnew(montecarlo_loop, :),NSCCdelays_usec2,AVGrates2,TOTALspikes2] = SCCfull_m({SpikeTrains_p, SpikeTrains_n},DELAYbinwidth,SCCdur);
%     [NSCCnew_debug,NSCCdelays_usec2_debug,AVGrates2_debug,TOTALspikes2_debug] = SCCfull_m({SpikeTrains_p, SpikeTrains_n},DELAYbinwidth,SCCdur);
    T(montecarlo_loop).new_sacc=toc;
end

figure(1);
clf;
histogram([T.old_sacc],100);
hold on;
histogram([T.new_sacc],100);
legend('old', 'new');
title(sprintf('old mean time --> %1.4f, new mean time --> %1.4f', mean([T.old_sacc]), mean([T.new_sacc])));


lw=1.5;
figure(2); clf;
plot(NSCCdelays_usec, mean(NSCCold,1), 'linewidth', lw)
hold on;
plot(NSCCdelays_usec, mean(NSCCnew,1), '--', 'linewidth', lw)
legend('old', 'new');