function [NSCC,NSCCdelays_usec,AVGrates,TOTALspikes] = SCCfull_m(SpikeTrains,DELAYbinwidth,Duration)
% File: ShufCrossCorr
% 27Sep2004: M Heinz - updated from Recruitment - NoiseCorr Analysis
%
% Created: 9/3/03 M. Heinz
% Calls: InnerSACmex.dll  (MEX-C file: InnerSACmex.c)
%
% Computes Normalized Shuffled Cross-Correlogram (NSCC) from a set of Spike Trains and Duration

NUMspikeREPS=cell(1,2);
for UNITind=1:2
    NUMspikeREPS{UNITind}=length(SpikeTrains{UNITind});
end

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
Kmax=cell(1,2);
NUMspikes=cell(1,2);
TOTALspikes=cell(1,2);
AVGrates=cell(1,2);
for UNITind=1:2
    NUMspikes{UNITind}=cellfun(@(x) numel(x),SpikeTrains{UNITind})';
    Kmax{UNITind}=max(NUMspikes{UNITind});
    TOTALspikes{UNITind}=sum(NUMspikes{UNITind});
    AVGrates{UNITind}=TOTALspikes{UNITind}/NUMspikeREPS{UNITind}/Duration;
end

%%%% Compute AVGrates
SpikeMAT=cell(1,2);

for UNITind=1:2
    SpikeMAT{UNITind}=NaN*ones(NUMspikeREPS{UNITind},Kmax{UNITind});
    for REPindREF=1:NUMspikeREPS{UNITind}
        SpikeMAT{UNITind}(REPindREF,1:length(SpikeTrains{UNITind}{REPindREF}))=SpikeTrains{UNITind}{REPindREF};
    end    
end

[SCC,~] = SCCfull(SpikeMAT{1}',NUMspikes{1},TOTALspikes{1},SpikeMAT{2}',NUMspikes{2},TOTALspikes{2},Duration, DELAYbinwidth);

NSCCdelays_usec=0:DELAYbinwidth:Duration;
%NSCCdelays=0:DELAYbinwidth:Duration/4;
NSCCdelays_usec=[-fliplr(NSCCdelays_usec(2:end)) NSCCdelays_usec]/1e-6;  % Convert into micro-seconds
% SCC=[fliplr(intsMEX(2:end)) intsMEX];  % Convert into micro-seconds

% SCC=hist(ints,NSCCdelays);
NSCC=SCC/(NUMspikeREPS{1}*NUMspikeREPS{2}*Duration*AVGrates{1}*AVGrates{2}*DELAYbinwidth);

return;