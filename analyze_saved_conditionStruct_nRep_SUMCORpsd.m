function analyze_saved_conditionStruct_nRep_SUMCORpsd(conditionStruct, stimFileNames, OUTDir)
% conditionStruct;
% stimFileNames;

uniq_stim=unique(stimFileNames);
uniq_SRtypes=unique([conditionStruct.SRtype]);
uniq_CFskHz=unique([conditionStruct.CF_kHz]);

uniq_nReps=unique([conditionStruct.nReps]);
uniq_windows=unique([conditionStruct.window]);
figOutDir=[OUTDir 'png_figs' filesep];
if ~isdir(figOutDir)
    mkdir(figOutDir);
end

for stimVar=1:length(uniq_stim)
    for srVar=1:length(uniq_SRtypes)
        curSR=uniq_SRtypes(srVar);
        for cfVar=1:length(uniq_CFskHz)
            curCF=uniq_CFskHz(cfVar);
            figFileName=sprintf('comp_sent%d_sr%d_cf%.0f', stimVar, curSR, curCF);
            clf;
            for nRepVar=1:length(uniq_nReps)
                curnRep=uniq_nReps(nRepVar);
                for windowVar=1:length(uniq_windows)
                    curWin=uniq_windows(windowVar);
                    curIND= [conditionStruct.stim]==stimVar & [conditionStruct.SRtype]==curSR & [conditionStruct.CF_kHz]==curCF & [conditionStruct.nReps]==curnRep & [conditionStruct.window]==curWin;
                    curIND=find(curIND==1);
                    if numel(curIND)~=1
                        error('Why multiple/zero match!');
                    else
                        subplot(length(uniq_nReps), length(uniq_windows), windowVar+(nRepVar-1)*length(uniq_nReps));
                        yyaxis left;
                        plot(conditionStruct(curIND).N_freqVEC, 20*log10(conditionStruct(curIND).N_PSDenv));
                        ylabel('power (dB)');
                        ylimLeft=ylim;
                        
                        yyaxis right;
                        plot(conditionStruct(curIND).uR_freqVEC, 20*log10(conditionStruct(curIND).uR_PSDenv));
                        xlim([0 300]);
                        xlabel('freq (Hz)');
                        ylabel('power (dB)');
                        ylimRight=ylim;
                        
                        if range(ylimLeft)>range(ylimRight)
                            yyaxis right;
                            ylim([mean(ylimRight)-.5*range(ylimLeft) mean(ylimRight)+.5*range(ylimLeft)]);
                        else 
                            yyaxis left;
                            ylim([mean(ylimLeft)-.5*range(ylimRight) mean(ylimLeft)+.5*range(ylimRight)]);
                        end
                        
                        title(strrep(sprintf('nRep%d / window%.3f', curnRep, curWin), '_', '|'));
                    end
                end
            end
            subplot(length(uniq_nReps), length(uniq_windows), 1);
            legend('spikes', 'uRate', 'location', 'southwest');
            set(gcf, 'units', 'normalized', 'position', [0 0 1 1]);
            set(gca, 'fontsize', 10);
            saveas(gcf, [figOutDir figFileName '.png'], 'png');
        end
    end
end