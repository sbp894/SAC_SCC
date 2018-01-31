%%
% -----------------------Results-------------------------
% 64 ms is pretty good in most cases for 25 repititions. SUMCOR PSDs
% estimated from spikes and meanrates closely match. 
% 
% Created by SP
% Testing
% Sumcor PSD convergence for L/M/H SR fibers (effect of duration and nrep)
% Naming convention:
% N* stands for neural
% uR* stand for meanrate

%%
clear;
close all;
clc;

addpath(['mexsource' filesep]);
OUTDir='effects_nrep_sumcorPSD_OUTDIR/';
if ~isdir(OUTDir)
   mkdir(OUTDir); 
end

OUTDirIter=[OUTDir 'singIter' filesep];
if ~isdir(OUTDirIter)
    mkdir(OUTDirIter);
end


%%
lw2=1.7;
lw1=1.2;
stimFileNames={'Stimuli/stimSetFluctuating/Stim0dB_SN_P.wav', 'Stimuli/stimSetStationary/Stim-6dB_SN_N.wav', 'Stimuli/stimSetStationary/Stim3dB_N_P.wav', 'Stimuli/stimSetFluctuating/Stim-3dB_N_N.wav'};
% stimFileNames={'Stimuli/stimSetFluctuating/Stim0dB_SN_P.wav'};
%%
all_CFs=[.5 1 2 4 8];
all_nReps=[15, 25, 1000];
all_SRtype=1:3;
all_windows=[.064 .128 .512]; % Seconds
all_stims=1:length(stimFileNames); % to use combVec directly
conditionVecs=combvec(all_CFs, all_nReps, all_SRtype, all_windows, all_stims)'; % All these should be row vectors to work, output is N by 5
conditionStruct=cell2struct(num2cell(conditionVecs), {'CF_kHz', 'nReps', 'SRtype', 'window', 'stim'}, 2);

srLetters='LMH';

AN.fs=100e3;
AN.Cohc=1;
AN.Cihc=1;
AN.species=2;
AN.noiseType=0;
AN.implnt=0;

params.DELAYbinwidth=50e-6;
params.SACCtrifiltWidth=5;
params.SCtrifiltWidth=1;
params.time_halfbandwidth=1;
params.num_seq = 2*params.time_halfbandwidth-1;
params.Fs_PSD=1/params.DELAYbinwidth;

parfor condVar=1:length(conditionStruct)
    fprintf('Start - %d\n', condVar);

    CF_Hz=conditionStruct(condVar).CF_kHz*1e3;
    nReps=conditionStruct(condVar).nReps;
    SRtype=conditionStruct(condVar).SRtype;
    window=conditionStruct(condVar).window;
    stimFname=stimFileNames(conditionStruct(condVar).stim); %#ok<*PFBNS>
    
    [stim, fsStim]=audioread(stimFname{1});
    windowCenter=window/2 + (length(stim)/fsStim-window)*rand(1);
    indStart=round(fsStim*(windowCenter-window/2));
    inEnd=round(fsStim*(windowCenter+window/2));
    stim_windowed=stim(indStart:inEnd);
    
    
    stimDur=length(stim_windowed)/fsStim;
    SCCdur=.8*stimDur;
    MAXdelay_ind=floor(SCCdur/params.DELAYbinwidth);
    
    vIHC_plus = model_IHC(stim_windowed.',CF_Hz,1,1/AN.fs,stimDur+0.001,AN.Cohc,AN.Cihc,AN.species);
    vIHC_minus = model_IHC(-stim_windowed.',CF_Hz,1,1/AN.fs,stimDur+0.001,AN.Cohc,AN.Cihc,AN.species);
    
    [meanrate_adap_plus,meanrate_unad_plus, ~] = model_Synapse(vIHC_plus,CF_Hz,1,1/AN.fs,SRtype,AN.noiseType ,AN.implnt);
    [meanrate_adap_minus,meanrate_unad_minus, ~] = model_Synapse(vIHC_minus,CF_Hz,1,1/AN.fs,SRtype,AN.noiseType ,AN.implnt);
    
    SpikeTrains_plus=get_sptimes(meanrate_unad_plus, AN.fs, nReps);
    SpikeTrains_minus=get_sptimes(meanrate_unad_minus, AN.fs, nReps);
    
    [NSAC_plus,~,~,~] = SAChalf_m(SpikeTrains_plus,params.DELAYbinwidth,SCCdur);
    [NSAC_minus,~,~,~] = SAChalf_m(SpikeTrains_minus,params.DELAYbinwidth,SCCdur);
    NSAC_avg=trifilt((NSAC_plus+NSAC_minus)/2, params.SACCtrifiltWidth);
    
    [NSCC_plus,Ndelays,~,~] = SCCfull_m({SpikeTrains_plus,SpikeTrains_minus},params.DELAYbinwidth,SCCdur);
    NSCC_minus=fliplr(NSCC_plus);
    NSCC_avg=trifilt((NSCC_plus+NSCC_minus)/2, params.SACCtrifiltWidth);
    
    ZEROind=dsearchn(Ndelays', 0);
    SACinds=(ZEROind-MAXdelay_ind:ZEROind+MAXdelay_ind);
    % SACdelays_usec=delays_usec(SACinds);
    % SAC limited-data Window CORRECTION
    TEMP=linspace(0,SCCdur/stimDur,ZEROind);
    TEMP1=1./(1-TEMP);
    WindowCORRECTION=[fliplr(TEMP1(2:end)) TEMP1];
    WindowCORRECTION(isinf(WindowCORRECTION))=1/eps;
    N_sumcor=trifilt((NSAC_avg+NSCC_avg)/2, params.SCtrifiltWidth).*WindowCORRECTION;
    
    [E, V]=dpss(numel(N_sumcor), params.time_halfbandwidth, params.num_seq);
    Eeven_N=E(:, 1:2:end);
    Veven_N=V(1:2:end);
    Nfft_psd_N=2^nextpow2(2.01*SCCdur*round(1/params.DELAYbinwidth));
    
    if params.num_seq>=3
        [Pxx_pmtm,F_pmtm] = pmtm(N_sumcor'-mean(N_sumcor),Eeven_N, Veven_N, Nfft_psd_N, params.Fs_PSD);
        inds_pmtm_to_fft=[1:(length(Pxx_pmtm)-1) (length(Pxx_pmtm)-1):-1:1];
        FFTtempN=sqrt(Pxx_pmtm(inds_pmtm_to_fft));
        N_freqVEC=F_pmtm(inds_pmtm_to_fft);
    elseif params.num_seq==1
        FFTtempN=mean(abs(fft((N_sumcor'-mean(N_sumcor)).*Eeven_N, Nfft_psd_N, 1))/length(N_sumcor), 2);
        N_freqVEC=(0:Nfft_psd_N-1)/Nfft_psd_N*params.Fs_PSD;
    end
    
    CF_indexN=dsearchn(N_freqVEC', CF_Hz); % use SACSCC_CF_Hz
    FFTadjN=zeros(size(FFTtempN));
    FFTadjN(1:CF_indexN)=FFTtempN(1:CF_indexN);
    FFTadjN((length(FFTtempN)-CF_indexN+1):end)=FFTtempN((length(FFTtempN)-CF_indexN+1):end); %keep negative freqs
    N_PSDenv=abs(FFTadjN);
%     conditionStruct(condVar).N_PSDenv=N_PSDenv;
%     conditionStruct(condVar).N_freqVEC=N_freqVEC;
    
%%    
    [uR_ACF_plus, ~]=xcorr(meanrate_unad_plus, 'unbiased');
    [uR_ACF_minus, ~]=xcorr(meanrate_unad_minus, 'unbiased');
    uR_ACF=(uR_ACF_plus+uR_ACF_minus)/2;
    [uR_CCF, delay_meanrate]=xcorr(meanrate_unad_minus, meanrate_unad_plus, 'unbiased');
    delay_meanrate=delay_meanrate/fsStim;
    uR_sumcor=resample((uR_ACF+uR_CCF)/2, params.Fs_PSD, fsStim);
    delay_meanrate=linspace(min(delay_meanrate), max(delay_meanrate), length(uR_sumcor));
    validInds= abs(delay_meanrate)<SCCdur;
    delay_meanrate=delay_meanrate(validInds);
    uR_sumcor=uR_sumcor(validInds);
    
    [E, V]=dpss(numel(uR_sumcor), params.time_halfbandwidth, params.num_seq);
    Eeven_uR=E(:, 1:2:end);
    Veven_uR=V(1:2:end);
    Nfft_psd_uR=2^nextpow2(2.01*SCCdur*round(1/params.DELAYbinwidth));
    if params.num_seq>=3
        [Pxx_pmtm,F_pmtm] = pmtm(uR_sumcor'-mean(uR_sumcor),Eeven_uR, Veven_uR, Nfft_psd_uR, params.Fs_PSD);
        inds_pmtm_to_fft=[1:(length(Pxx_pmtm)-1) (length(Pxx_pmtm)-1):-1:1];
        FFTtempuR=sqrt(Pxx_pmtm(inds_pmtm_to_fft));
        uR_freqVEC=F_pmtm(inds_pmtm_to_fft);
    elseif params.num_seq==1
        FFTtempuR=mean(abs(fft((uR_sumcor'-mean(uR_sumcor)).*Eeven_uR, Nfft_psd_uR, 1))/length(uR_sumcor), 2);
        uR_freqVEC=(0:Nfft_psd_uR-1)/Nfft_psd_uR*params.Fs_PSD;
    end
    CF_indexuR=dsearchn(uR_freqVEC', CF_Hz); % use SACSCC_CF_Hz
    uR_FFTadj=zeros(size(FFTtempuR));
    uR_FFTadj(1:CF_indexuR)=FFTtempuR(1:CF_indexuR);
    uR_FFTadj((length(FFTtempuR)-CF_indexuR+1):end)=FFTtempuR((length(FFTtempuR)-CF_indexuR+1):end); %keep negative freqs
    uR_PSDenv=abs(uR_FFTadj);
%     conditionStruct(condVar).uR_PSDenv=uR_PSDenv;
%     conditionStruct(condVar).uR_freqVEC=uR_freqVEC;
    
    dummyFunforSaving([OUTDirIter 'conditionData_iter_' num2str(condVar) '.mat'], uR_freqVEC, uR_PSDenv, N_freqVEC, N_PSDenv, CF_Hz, nReps, SRtype, window, stimFname);
end

% save([OUTDir 'conditionData.mat'], 'conditionStruct', 'stimFileNames', 'AN', 'params');
% analyze_saved_conditionStruct_nRep_SUMCORpsd(conditionStruct, stimFileNames, OUTDir);
run_analyze_saved_conditionStruct_nRep_SUMCORpsd;