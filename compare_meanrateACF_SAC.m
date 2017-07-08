%%
% Created by SP
% Testing 
% (1) Whether meanrate ACF is same as SAC
% (2) The effects of nReps on when SAC converges

%%
clear;
close all;
clc;

addpath('mexsource\');


%%
lw2=1.7;
lw1=1.2;

%%
all_nReps=[10, 15, 20, 100];
all_SRtype=1:3;
srLetters='LMH';

CF_Hz=1e3;
ANmodel_Fs_Hz=100e3;
Cohc=1;
Cihc=1;
species=2;
noiseType=0;
implnt=0;
DELAYbinwidth=50e-6;

[stim, fs]=audioread('dan_sent1.wav');
stimDur=length(stim)/fs;
Duration=.8*stimDur;



corr_ACFs=nan(length(all_SRtype), length(all_nReps));

for srVar=1:length(all_SRtype)
    SRtype=all_SRtype(srVar);
    
    for nRepVar=1:length(all_nReps)
        
        nReps=all_nReps(nRepVar);
        
        
        %%
        
        
        vIHC = model_IHC(stim.',CF_Hz,1,1/ANmodel_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
        [meanrate_adap,meanrate_unad, ~] = model_Synapse(vIHC,CF_Hz,1,1/ANmodel_Fs_Hz,SRtype,noiseType ,implnt);
        
        SpikeTrains=get_sptimes(meanrate_unad, ANmodel_Fs_Hz, nReps);
        [NSAC,NSACdelays,AVGrates,TOTALspikes] = SAChalf_m(SpikeTrains,DELAYbinwidth,Duration);
        WindowCorrection=1./(stimDur-abs(NSACdelays));
        NSAC=NSAC.*WindowCorrection;
        NSAC=trifilt(NSAC, 10);
        
        [ACF_meanrate, delay_meanrate]=xcorr(meanrate_adap, 'unbiased');
        delay_meanrate=delay_meanrate/ANmodel_Fs_Hz;
        ACF_meanrate=ACF_meanrate*max(NSAC)/max(ACF_meanrate);
        
        
        %%
        if rem(ANmodel_Fs_Hz*DELAYbinwidth,1)==0
            ACF_meanrate_rs=resample(ACF_meanrate, 1, ANmodel_Fs_Hz*DELAYbinwidth);
            zeroInd=(length(ACF_meanrate_rs)+1)/2;
            indHalf=(length(NSACdelays)-1)/2;
            ACF_meanrate_rs=ACF_meanrate_rs(zeroInd-indHalf:zeroInd+indHalf);
            corr_ACFs(srVar, nRepVar)=corr(ACF_meanrate_rs', NSAC');
            %%
            figure(101);
            subplot(length(all_SRtype), length(all_nReps), (srVar-1)*length(all_nReps)+nRepVar);
            plot(NSACdelays, ACF_meanrate_rs);
            hold on;
            plot(NSACdelays, NSAC);
            %             xlim([-.05 .05]);
            title(sprintf('Corr=%.2f::%sSR: nRep=%i',corr_ACFs(srVar, nRepVar), srLetters(SRtype), nReps));
            
        else
            error('binwidth*AN_fs should be integer\n');
        end
        
        %%
%         figure(srVar);
%         subplot(length(all_nReps),1,nRepVar);
%         plot(delay_meanrate, ACF_meanrate, 'linewidth', lw2);
%         hold on;
%         plot(NSACdelays, NSAC,'linewidth', lw1);
%         xlim([-.05 .05]);
%         title(sprintf('%sSR: nRep=%i',srLetters(SRtype), nReps));
        
    end
    xlabel('delay: seconds');
end