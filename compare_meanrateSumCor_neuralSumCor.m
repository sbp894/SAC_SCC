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
[stim, fs]=audioread('Stim_S_P.wav');
% [stim, fs]=audioread('Stim-6dB_N_P.wav');
spl2use=65;
stim=stim*(20e-6)*10^(spl2use/20)/rms(stim);
stimDur=length(stim)/fs;
Duration=.4*stimDur;

%%
h_pMod=234;
h_SC=746;

figure(h_pMod);
set(h_pMod, 'visible', 'off');
figure(h_SC);
set(h_SC, 'visible', 'off');

lw2=1.5;
lw1=1.2;
lw=2;

markers='dospv';


%%
all_nReps=[10 20 50];
all_SRtype=1:3;
srLetters='LMH';

%%
CF_Hz=.5e3;
AN_Fs_Hz=100e3;
Cohc=1;
Cihc=1;
species=2;
noiseType=0;
implnt=0;
DELAYbinwidth=50e-6;
SC_Fs_Hz=1/DELAYbinwidth;

%%
if AN_Fs_Hz~=fs
   error('sampling rates are not equal for model and sitmulus'); 
end

%%
d=fdesign.lowpass('Fp,Fst,Ap,Ast',1.8*CF_Hz,2.2*CF_Hz,.1, 60, SC_Fs_Hz);
Hd=design(d, 'cheby2');
[b,a]=sos2tf(Hd.sosMatrix);

corr_SACs=nan(length(all_SRtype), length(all_nReps));
corr_SCs=nan(length(all_SRtype), length(all_nReps));

plt=0;
ModFreqs=2.^(0:6);
qFactor=ones(size(ModFreqs));
strModFreqs=textscan(num2str(ModFreqs), '%s');
strModFreqs=strModFreqs{1};

legstr=cell(length(all_nReps)+2, 1);
legstr{1}='time-adap';
legstr{2}='time-unad';
for nRepVar=1:length(all_nReps)
    legstr{nRepVar+2}=sprintf('nRep=%i', all_nReps(nRepVar));
end

for srVar=1:length(all_SRtype)
    SRtype=all_SRtype(srVar);
    
    for nRepVar=1:length(all_nReps)
        
        nReps=all_nReps(nRepVar);
        
        
        %% For positive polarity
        vIHCplus= model_IHC(stim.',CF_Hz,1,1/AN_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
        [meanrate_adap_plus,meanrate_unad_plus, ~]= model_Synapse(vIHCplus,CF_Hz,1,1/AN_Fs_Hz,SRtype,noiseType ,implnt);
        SpikeTrainsplus=get_sptimes(meanrate_unad_plus, AN_Fs_Hz, nReps);
        
        
        [n_NSAC_plus,NSACdelays,~,~] = SAChalf_m(SpikeTrainsplus,DELAYbinwidth,Duration);
        n_WindowCorrection=stimDur./(stimDur-abs(NSACdelays));
        n_NSAC_plus=n_NSAC_plus.*n_WindowCorrection;
        n_NSAC_plus=trifilt(n_NSAC_plus, 10);
        
        %% For negative polarity
        vIHCminus= model_IHC(stim.',CF_Hz,1,1/AN_Fs_Hz,stimDur+0.05,Cohc,Cihc,species);
        [meanrate_adap_minus,meanrate_unad_minus, ~]= model_Synapse(vIHCminus,CF_Hz,1,1/AN_Fs_Hz,SRtype,noiseType ,implnt);
        SpikeTrainsminus=get_sptimes(meanrate_unad_minus, AN_Fs_Hz, nReps);
        
        
        [n_NSAC_minus,NSACdelays,~,~] = SAChalf_m(SpikeTrainsminus,DELAYbinwidth,Duration);
        n_WindowCorrection=stimDur./(stimDur-abs(NSACdelays));
        n_NSAC_minus=n_NSAC_minus.*n_WindowCorrection;
        n_NSAC_minus=trifilt(n_NSAC_minus, 10);
        
        n_NSAC=(n_NSAC_plus+n_NSAC_minus)/2;
        
        %%
        [n_NSCC,NSCCdelays,AVGrates,TOTALspikes] = SCCfull_m({SpikeTrainsplus, SpikeTrainsminus},DELAYbinwidth,Duration);
        n_WindowCorrection=stimDur./(stimDur-abs(NSCCdelays));
        n_NSCC=n_NSCC.*n_WindowCorrection;
        n_NSCC=trifilt(n_NSCC, 10);
        
        n_SUMCOR=(n_NSAC+n_NSCC)/2;
        n_SUMCOR=filtfilt(b, a, n_SUMCOR);
        
        [t_NSAC_plus, delay_meanrate]=xcorr(meanrate_adap_plus, 'coeff');
        t_WindowCorrection=length(meanrate_adap_plus)./(length(meanrate_adap_plus)-abs(delay_meanrate));
        t_NSAC_plus=t_NSAC_plus.*t_WindowCorrection;
        [t_NSAC_minus, ~]=xcorr(meanrate_adap_minus, 'coeff');
        t_NSAC_minus=t_NSAC_minus.*t_WindowCorrection;
        t_NSAC=(t_NSAC_plus+t_NSAC_minus)/2;
        
        
        [CCF_meanrate, ~]=xcorr(meanrate_adap_plus, meanrate_adap_minus, 'coeff');
        CCF_meanrate=CCF_meanrate.*t_WindowCorrection;
        t_ad_SUMCOR=(t_NSAC+CCF_meanrate)/2;
        
        
        [t_ua_NSAC_plus, ~]=xcorr(meanrate_unad_plus, 'coeff');
        [t_ua_NSAC_minus, ~]=xcorr(meanrate_unad_minus, 'coeff');
        t_ua_NSAC=(t_ua_NSAC_plus+t_ua_NSAC_minus)/2;
        t_ua_NSAC=t_ua_NSAC.*t_WindowCorrection;
        
        [CCF_ua_meanrate, ~]=xcorr(meanrate_unad_plus, meanrate_unad_minus, 'coeff');
        CCF_ua_meanrate=CCF_ua_meanrate.*t_WindowCorrection;
        t_ua_SUMCOR=(t_ua_NSAC+CCF_ua_meanrate)/2;
        delay_meanrate=delay_meanrate/AN_Fs_Hz;
        
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
            t_ua_SUMCOR_rs=resample(t_ua_SUMCOR, 1, AN_Fs_Hz/SC_Fs_Hz);
            
            zeroInd=round((length(t_NSAC_rs)+1)/2);
            indHalf=(length(NSACdelays)-1)/2;
            t_NSAC_rs=t_NSAC_rs(zeroInd-indHalf:zeroInd+indHalf);
            t_ad_SUMCOR_rs=t_ad_SUMCOR_rs(zeroInd-indHalf:zeroInd+indHalf);
            t_ua_SUMCOR_rs=t_ua_SUMCOR_rs(zeroInd-indHalf:zeroInd+indHalf);
            
            %% Normalize to same numbers.
            t_NSAC_rs=t_NSAC_rs*max(n_NSAC)/max(t_NSAC_rs);
            t_ad_SUMCOR_rs=t_ad_SUMCOR_rs*max(n_SUMCOR)/max(t_ad_SUMCOR_rs);
            t_ua_SUMCOR_rs=t_ua_SUMCOR_rs*max(n_SUMCOR)/max(t_ua_SUMCOR_rs);
            
            corr_SACs(srVar, nRepVar)=corr(t_NSAC_rs', n_NSAC_plus');
            corr_SCs(srVar, nRepVar)=corr(t_ad_SUMCOR_rs', n_SUMCOR');
            
            %%
            NFFT=2^nextpow2(length(n_SUMCOR));
            PSDfreqVEC_Hz=linspace(-SC_Fs_Hz/2, SC_Fs_Hz/2, NFFT)';
            binWeights = modFbank(PSDfreqVEC_Hz,ModFreqs, qFactor, plt);
            
            n_SC_PSD=abs(fftshift(fft(n_SUMCOR-mean(n_SUMCOR), NFFT)))/length(n_SUMCOR);
            t_ad_SC_PSD=abs(fftshift(fft(t_ad_SUMCOR_rs-mean(t_ad_SUMCOR_rs), NFFT)))/length(t_ad_SUMCOR_rs);
            t_ua_SC_PSD=abs(fftshift(fft(t_ua_SUMCOR_rs-mean(t_ua_SUMCOR_rs), NFFT)))/length(t_ua_SUMCOR_rs);
            
            n_PowerMod= nansum(binWeights.* repmat(n_SC_PSD, length(ModFreqs), 1),2);
            t_ad_PowerMod= nansum(binWeights.* repmat(t_ad_SC_PSD, length(ModFreqs), 1),2);
            t_ua_PowerMod= nansum(binWeights.* repmat(t_ua_SC_PSD, length(ModFreqs), 1),2);
            
            
            %%
            %             figure(101);
            %             subplot(length(all_SRtype), length(all_nReps), (srVar-1)*length(all_nReps)+nRepVar);
            %             plot(NSACdelays, ACF_meanrate_rs);
            %             hold on;
            %             plot(NSACdelays, NSAC);
            %             %             xlim([-.05 .05]);
            %             title(sprintf('Corr=%.2f::%sSR: nRep=%i',corr_ACFs(srVar, nRepVar), srLetters(SRtype), nReps));
            %
            
            %%
            figure(h_SC);
            subplot(length(all_SRtype), length(all_nReps), (srVar-1)*length(all_nReps)+nRepVar);
            hold on;
            plot(NSACdelays, t_ad_SUMCOR_rs, 'k', 'linewidth', lw1);
            plot(NSACdelays, t_ua_SUMCOR_rs, '--m', 'linewidth', lw1);
            plot(NSACdelays, n_SUMCOR, 'linewidth', lw1);
            title(sprintf('Corr=%.2f::%sSR: nRep=%i',corr_SCs(srVar, nRepVar), srLetters(SRtype), nReps));
            %             xlim([-.05 .05]);
            
            %%
            figure(h_pMod);
            subplot(length(all_SRtype), 1, srVar);
            if nRepVar==1
                hold on;
                plot(t_ad_PowerMod, '-*k', 'linewidth', lw);
                plot(t_ua_PowerMod, '--^m', 'linewidth', lw);
            end
            plot(n_PowerMod, ['--' markers(nRepVar)], 'linewidth', lw);
            axis tight;
            set(gca, 'xticklabel',  strModFreqs);
            %             xlim([-.05 .05]);
            
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
    figure(h_pMod);
    title(sprintf('%sSR',srLetters(SRtype)));
end
figure(h_SC);
legend('time- adapted', 'time- unadapted', 'neural');
xlabel('delay: seconds');

figure(h_pMod);
legend(legstr);
xlabel('modulation freqs');

set(h_pMod, 'visible', 'on');
set(h_SC, 'visible', 'on');
