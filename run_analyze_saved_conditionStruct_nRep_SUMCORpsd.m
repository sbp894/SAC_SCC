clear;
close all;
clc;


%% ----------------------Block Starts-----------------------------------
% Need to update! Copied it directly from the effects_nrep_sumcorPSD

addpath(['mexsource' filesep]);
OUTDir='effects_nrep_sumcorPSD_OUTDIR/';
if ~isdir(OUTDir)
   mkdir(OUTDir); 
end

OUTDirIter=[OUTDir 'singIter' filesep];
if ~isdir(OUTDirIter)
    mkdir(OUTDirIter);
end


%
lw2=1.7;
lw1=1.2;
stimFileNames={'Stimuli/stimSetFluctuating/Stim0dB_SN_P.wav', 'Stimuli/stimSetStationary/Stim-6dB_SN_N.wav', 'Stimuli/stimSetStationary/Stim3dB_N_P.wav', 'Stimuli/stimSetFluctuating/Stim-3dB_N_N.wav'};
% stimFileNames={'Stimuli/stimSetFluctuating/Stim0dB_SN_P.wav'};
%
all_CFs=[.5 1 2 4 8];
all_nReps=[15, 25, 1000];
all_SRtype=1:3;
all_windows=[.064 .128 .512]; % Seconds
all_stims=1:length(stimFileNames); % to use combVec directly
conditionVecs=combvec(all_CFs, all_nReps, all_SRtype, all_windows, all_stims)'; % All these should be row vectors to work, output is N by 5
cSStruct=cell2struct(num2cell(conditionVecs), {'CF_kHz', 'nReps', 'SRtype', 'window', 'stim'}, 2);
% ----------------------Block Ends-----------------------------------

%%

DATADir= 'effects_nrep_sumcorPSD_OUTDIR/singIter/';
condStruct=load('effects_nrep_sumcorPSD_OUTDIR/conditionStruct.mat');
allfiles=dir([DATADir '*.mat']);

for fileVar=1:length(allfiles)
   data=load([DATADir allfiles(fileVar).name]);
   structFields=fieldnames(data); 
   for fieldVar=1:length(structFields)
       eval(sprintf('conditionStruct(%d).%s=data.%s;', fileVar, structFields{fieldVar}, structFields{fieldVar}));
   end   
   %% Should be unnecessary if rerun and save data with updated effect_nrep_sumcorPSD 
   % ----------------------Block Starts-------------------------------------------
   iterNum=sscanf(allfiles(fileVar).name, 'conditionData_iter_%d.mat');
   structFields=fieldnames(cSStruct(iterNum)); 
   for fieldVar=1:length(structFields)
       eval(sprintf('conditionStruct(%d).%s=cSStruct(iterNum).%s;', fileVar, structFields{fieldVar}, structFields{fieldVar}));
   end 
   % ------------------------Block Ends-------------------------------------------
end

analyze_saved_conditionStruct_nRep_SUMCORpsd(conditionStruct, stimFileNames, OUTDir);