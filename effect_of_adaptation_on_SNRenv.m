%%
% Created by SP
% To see the effects of adaptation on SUMCOR PSD. Adaptation acts as a low
% pass filter. So adaptation should results in higher power in lower
% modulation frequencies. We are getting higher power using neural sEPSM in
% 1 and 2 Hz modulation frequency bands.

%%
clear;
close all;
clc;

addpath('mexsource\');

%%
all_nReps=5;
all_SRtype=1:3;
SNR=-6;
All_stims={{['Stimuli\stimSetStationary\Stim' num2str(SNR) 'dB_SN_P.wav'], ['Stimuli\stimSetStationary\Stim' num2str(SNR) 'dB_N_P.wav']}, ...
    {['Stimuli\stimSetFluctuating\Stim' num2str(SNR) 'dB_SN_P.wav'], ['Stimuli\stimSetFluctuating\Stim' num2str(SNR) 'dB_N_P.wav']}};
for stimVar=1:length(All_stims)
    stimNames=All_stims{stimVar};
    srLetters='LMH';
    
    %%
    AN.CF_Hz=1e3;
    AN.AN_Fs_Hz=100e3;
    AN.Cohc=1;
    AN.Cihc=1;
    AN.species=2;
    AN.noiseType=0;
    AN.implnt=0;
    
    %%
    params.nPointsSmooth=10;
    params.ModFreqs=2.^(0:6);
    params.strModFreqs=textscan(num2str(params.ModFreqs),'%s');
    params.strModFreqs=params.strModFreqs{1};
    params.spl2use=65;
    params.maxDelayFrac=.6;
    params.DELAYbinwidth=50e-6;
    params.plotYes=0;
    params.useNeural=0;
    params.useNeural=0;
    
    %%
    hFig.lw=2;
    hFig.h_pMod=234;
    hFig.h_SC=746;
    hFig.SNRenv=1;
    hFig.linestyles={'-', '--', '-.'};
    figure(stimVar);
    co=get(gca, 'colororder');
    
    %%
    n_SNRenv=nan(length(all_SRtype), length(all_nReps), length(params.ModFreqs));
    t_ad_SNRenv=nan(length(all_SRtype), length(all_nReps), length(params.ModFreqs));
    t_ua_SNRenv=nan(length(all_SRtype), length(all_nReps), length(params.ModFreqs));
    
    if params.useNeural
        legstr=cell(length(all_nReps)+2,1);
        legstr2=cell(2*(length(all_nReps)+2),1);
    else
        legstr=cell(2,1);
        legstr2=cell(4,1);
    end
    legstr(1:2)={'tAdapt', 'tUnad'};
    legstr2(1:4)={'SN-tAdapt', 'SN-tUnad','N-tAdapt', 'N-tUnad'};
    if params.useNeural
        for nRepVar=1:length(all_nReps)
            legstr{nRepVar+2}=sprintf('nRep=%i', all_nReps(nRepVar));
            legstr2{2*(nRepVar+2)-1}=sprintf('SN-nRep=%i', all_nReps(nRepVar));
            legstr2{2*(nRepVar+2)}=sprintf('N-nRep=%i', all_nReps(nRepVar));
        end
    end
    %% Main loop
    for srVar=1:length(all_SRtype)
        params.SRtype=all_SRtype(srVar);
        
        for nRepVar=1:length(all_nReps)
            
            params.nReps=all_nReps(nRepVar);
            
            % SN
            params.StimName=stimNames{1};
            [n_PowerMod_SN, t_ad_PowerMod_SN, t_ua_PowerMod_SN]=fun_getpMOD(params, AN, hFig);
            subplot(length(all_SRtype),2,2*srVar-1);
            hold on;
            if nRepVar==1
                plot(t_ad_PowerMod_SN, [ hFig.linestyles{1} 's'], 'color', co(1,:), 'linewidth', hFig.lw);
                plot(t_ua_PowerMod_SN, [ hFig.linestyles{1} 's'], 'color', co(2,:), 'linewidth', hFig.lw);
            end
            if params.useNeural
                plot(n_PowerMod_SN, [ hFig.linestyles{1} 'o'], 'color', co(2+nRepVar,:), 'linewidth', hFig.lw);
            end
            % N
            params.StimName=stimNames{2};
            [n_PowerMod_N, t_ad_PowerMod_N, t_ua_PowerMod_N]=fun_getpMOD(params, AN, hFig);
            subplot(length(all_SRtype),2,2*srVar-1);
            hold on;
            if nRepVar==1
                plot(t_ad_PowerMod_N, [ hFig.linestyles{2} 'd'], 'color', co(1,:), 'linewidth', hFig.lw);
                plot(t_ua_PowerMod_N, [ hFig.linestyles{2} 'd'], 'color', co(2,:), 'linewidth', hFig.lw);
            end
            if params.useNeural
                plot(n_PowerMod_N, [ hFig.linestyles{2} 'o'], 'color', co(2+nRepVar,:), 'linewidth', hFig.lw);
            end
            
            set(gca, 'XTick', 1:7, 'xticklabel',  params.strModFreqs);
            
            %% SNRenv calculations
            if params.useNeural
                nSe=(n_PowerMod_SN-n_PowerMod_N)./n_PowerMod_N;
                nSe(nSe<0)=0;
                n_SNRenv(srVar, nRepVar, :)=nSe;
            end
            
            t_ad_Se=(t_ad_PowerMod_SN-t_ad_PowerMod_N)./t_ad_PowerMod_N;
            t_ad_Se(t_ad_Se<0)=0;
            t_ad_SNRenv(srVar, nRepVar, :)=t_ad_Se;
            
            t_ua_Se=(t_ua_PowerMod_SN-t_ua_PowerMod_N)./t_ua_PowerMod_N;
            t_ua_Se(t_ua_Se<0)=0;
            t_ua_SNRenv(srVar, nRepVar, :)=t_ua_Se;
            
            %% Plot
            subplot(length(all_SRtype),2,2*srVar);
            hold on;
            if nRepVar==1
                plot(t_ad_Se, '-s', 'linewidth', hFig.lw);
                plot(t_ua_Se, '-s', 'linewidth', hFig.lw);
            end
            if params.useNeural
                plot(nSe, '-o', 'linewidth', hFig.lw);
            end
        end
        title(sprintf('%sSR, SNRenv=%.2f (t-Ad) and %.2f (t-unad)', srLetters(params.SRtype), nansum(t_ad_Se.^2), nansum(t_ua_Se.^2)));
        set(gca, 'XTick', 1:7, 'xticklabel',  params.strModFreqs);
    end
    
    legend(legstr, 'location', 'best');
    subplot(length(all_SRtype),2,2*srVar-1);
    legend(legstr2, 'location', 'best');
end