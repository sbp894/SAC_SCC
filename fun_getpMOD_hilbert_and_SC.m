function [n_PowerMod, t_ad_PowerMod, t_hilbert_PowerMod]=fun_getpMOD_hilbert_and_SC(params, AN, hFig)
n_PowerMod=[];

%%
[stim, fs]=audioread(params.StimName);
% [stim, fs]=audioread('Stim-6dB_N_P.wav');
spl2use=params.spl2use;
stim=stim*(20e-6)*10^(spl2use/20)/rms(stim);
stimDur=length(stim)/fs;
Duration=params.maxDelayFrac*stimDur;
plotYes=params.plotYes;
ModFreqs=params.ModFreqs;
DELAYbinwidth=params.DELAYbinwidth;
SC_Fs_Hz=1/DELAYbinwidth;

%%
if plotYes
    h_pMod=hFig.h_pMod;
    h_SC=hFig.h_SC;
    
    figure(h_pMod);
    set(h_pMod, 'visible', 'off');
    figure(h_SC);
    set(h_SC, 'visible', 'off');
    
    %     lw2=1.5;
    lw1=1.2;
    lw=2;
    
    markers='d';
end

srLetters='LMH';

%%
CF_Hz=AN.CF_Hz;
AN_Fs_Hz=AN.AN_Fs_Hz;
Cohc=AN.Cohc;
Cihc=AN.Cihc;
species=AN.species;
noiseType=AN.noiseType;
implnt=AN.implnt;



%%
if AN_Fs_Hz~=fs
    error('sampling rates are not equal for model and sitmulus');
end

%%
d=fdesign.lowpass('Fp,Fst,Ap,Ast',min(1.9*CF_Hz,SC_Fs_Hz/2.01),min(2.1*CF_Hz, SC_Fs_Hz/2),.1, 60, SC_Fs_Hz);
Hd=design(d, 'cheby2');
[b,a]=sos2tf(Hd.sosMatrix);

plt=0;
qFactor=ones(size(ModFreqs));
strModFreqs=textscan(num2str(ModFreqs), '%s');
strModFreqs=strModFreqs{1};


SRtype=params.SRtype;
nReps=params.nReps;
legstr={'time-adap','time-unad',sprintf('nRep=%i', nReps)};

%% For positive polarity

vIHCplus= model_IHC(stim.',CF_Hz,1,1/AN_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
[meanrate_adap_plus,meanrate_unad_plus, ~]= model_Synapse(vIHCplus,CF_Hz,1,1/AN_Fs_Hz,SRtype,noiseType ,implnt);
env_meanrate_adap_plus=abs(hilbert(meanrate_adap_plus));

if params.useNeural
    SpikeTrainsplus=get_sptimes(meanrate_unad_plus, AN_Fs_Hz, nReps);
    [n_NSAC_plus,NSACdelays,~,~] = SAChalf_m(SpikeTrainsplus,DELAYbinwidth,Duration);
    n_WindowCorrection=stimDur./(stimDur-abs(NSACdelays));
    n_NSAC_plus=n_NSAC_plus.*n_WindowCorrection;
    n_NSAC_plus=trifilt(n_NSAC_plus, params.nPointsSmooth);
end

%% For negative polarity
vIHCminus= model_IHC(stim.',CF_Hz,1,1/AN_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
[meanrate_adap_minus,meanrate_unad_minus, ~]= model_Synapse(vIHCminus,CF_Hz,1,1/AN_Fs_Hz,SRtype,noiseType ,implnt);


if params.useNeural
    SpikeTrainsminus=get_sptimes(meanrate_unad_minus, AN_Fs_Hz, nReps);
        
    [n_NSAC_minus,NSACdelays,~,~] = SAChalf_m(SpikeTrainsminus,DELAYbinwidth,Duration);
    n_WindowCorrection=stimDur./(stimDur-abs(NSACdelays));
    n_NSAC_minus=n_NSAC_minus.*n_WindowCorrection;
    n_NSAC_minus=trifilt(n_NSAC_minus, params.nPointsSmooth);
    n_NSAC=(n_NSAC_plus+n_NSAC_minus)/2;
    
    %%
    [n_NSCC,NSCCdelays,~,~] = SCCfull_m({SpikeTrainsplus, SpikeTrainsminus},DELAYbinwidth,Duration);
    n_WindowCorrection=stimDur./(stimDur-abs(NSCCdelays));
    n_NSCC=n_NSCC.*n_WindowCorrection;
    n_NSCC=trifilt(n_NSCC, params.nPointsSmooth);
    
    n_SUMCOR=(n_NSAC+n_NSCC)/2;
    n_SUMCOR=filtfilt(b, a, n_SUMCOR);
end

[t_NSAC_plus, delay_meanrate]=xcorr(meanrate_adap_plus, 'coeff');
t_WindowCorrection=length(meanrate_adap_plus)./(length(meanrate_adap_plus)-abs(delay_meanrate));
t_NSAC_plus=t_NSAC_plus.*t_WindowCorrection;
[t_NSAC_minus, ~]=xcorr(meanrate_adap_minus, 'coeff');
t_NSAC_minus=t_NSAC_minus.*t_WindowCorrection;
t_NSAC=(t_NSAC_plus+t_NSAC_minus)/2;


[CCF_meanrate, ~]=xcorr(meanrate_adap_plus, meanrate_adap_minus, 'coeff');
CCF_meanrate=CCF_meanrate.*t_WindowCorrection;
t_ad_SUMCOR=(t_NSAC+CCF_meanrate)/2;

% delay_meanrate=delay_meanrate/AN_Fs_Hz;

%%
%         figure(978);
%         clf;
%         plot(delay_meanrate, t_ad_SUMCOR);
%         hold on;
%         plot(delay_meanrate, t_ua_SUMCOR);
%         legend('adapted', 'unadapted');
%         xlabel('delay (sec)');
%         title('SUMCOR comparison for meanrates');
%%
if rem(AN_Fs_Hz/SC_Fs_Hz,1)==0
    t_NSAC_rs=resample(t_NSAC, 1, AN_Fs_Hz/SC_Fs_Hz);
    t_ad_SUMCOR_rs=resample(t_ad_SUMCOR, 1, AN_Fs_Hz/SC_Fs_Hz);
    t_env_hilbert_rs=resample(env_meanrate_adap_plus, 1, AN_Fs_Hz/SC_Fs_Hz);
    
    zeroInd=round((length(t_NSAC_rs)+1)/2);
    indHalf=floor(Duration/DELAYbinwidth);
%     t_NSAC_rs=t_NSAC_rs(zeroInd-indHalf:zeroInd+indHalf);
    t_ad_SUMCOR_rs=t_ad_SUMCOR_rs(zeroInd-indHalf:zeroInd+indHalf);
%     t_ua_SUMCOR_rs=t_ua_SUMCOR_rs(zeroInd-indHalf:zeroInd+indHalf);
    
    %% Normalize to same numbers.
%     t_NSAC_rs=t_NSAC_rs*max(n_NSAC)/max(t_NSAC_rs);
%     t_ad_SUMCOR_rs=t_ad_SUMCOR_rs*max(n_SUMCOR)/max(t_ad_SUMCOR_rs);
%     t_ua_SUMCOR_rs=t_ua_SUMCOR_rs*max(n_SUMCOR)/max(t_ua_SUMCOR_rs);
    
%     corr_SACs=corr(t_NSAC_rs', n_NSAC_plus');
    %%
    NFFT=2^nextpow2(length(t_ad_SUMCOR_rs));
    PSDfreqVEC_Hz=linspace(-SC_Fs_Hz/2, SC_Fs_Hz/2, NFFT)';
    binWeights = modFbank(PSDfreqVEC_Hz,ModFreqs, qFactor, plt);
    
    if params.useNeural
        corr_SCs=corr(t_ad_SUMCOR_rs', n_SUMCOR');
        n_SC_PSD=abs(fftshift(fft(n_SUMCOR-mean(n_SUMCOR), NFFT)))/length(n_SUMCOR);
        n_PowerMod= nansum(binWeights.* repmat(n_SC_PSD, length(ModFreqs), 1),2);
    end
    
    t_ad_SC_PSD=abs(fftshift(fft(t_ad_SUMCOR_rs-mean(t_ad_SUMCOR_rs), NFFT)))/length(t_ad_SUMCOR_rs);
    t_env_hilbert_PSD=abs(fftshift(fft(t_env_hilbert_rs-mean(t_env_hilbert_rs), NFFT)))/length(t_env_hilbert_rs);
    
    t_ad_PowerMod= nansum(binWeights.* repmat(t_ad_SC_PSD, length(ModFreqs), 1),2);
    t_hilbert_PowerMod= nansum(binWeights.* repmat(t_env_hilbert_PSD, length(ModFreqs), 1),2);
    
    %%
    if plotYes
        figure(h_SC);
        hold on;
        plot(NSACdelays, t_ad_SUMCOR_rs, 'k', 'linewidth', lw1);
        plot(NSACdelays, t_env_hilbert_rs, '--m', 'linewidth', lw1);
        plot(NSACdelays, n_SUMCOR, 'linewidth', lw1);
        title(sprintf('Corr=%.2f::%sSR: nRep=%i',corr_SCs, srLetters(SRtype), nReps));
        %             xlim([-.05 .05]);
        
        %%
        figure(h_pMod);
        
            hold on;
            plot(t_ad_PowerMod, '-*k', 'linewidth', lw);
            plot(t_hilbert_PowerMod, '--^m', 'linewidth', lw);
        plot(n_PowerMod, ['--' markers], 'linewidth', lw);
        axis tight;
        set(gca, 'xticklabel',  strModFreqs);
        %             xlim([-.05 .05]);
    end
    
else
    error('binwidth*AN_fs should be integer\n');
end

if plotYes
    figure(h_pMod);
    title(sprintf('%sSR',srLetters(SRtype)));
end

if plotYes
    figure(h_SC);
    legend('time- adapted', 'time- unadapted', 'neural');
    xlabel('delay: seconds');
    
    figure(h_pMod);
    legend(legstr);
    xlabel('modulation freqs');
    
    set(h_pMod, 'visible', 'on');
    set(h_SC, 'visible', 'on');
end
