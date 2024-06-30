function [coherence,freqs,binSteps] = coherence_AS(LFP1,LFP2,freqRange,timeWindow,stepWindowSize,plotIT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the time-frequency coherence between two signals.
% Thus, it breaks the signal into 2 second windows that are %50 overlapping
% and calculates the coherence during that time. The step windows are
% Hamming windows in order to deal with edge effects from the FFT. Also, it
% handles NaN values. If a NaN is detected within your signal, it analyzes
% the chucks of data surrounding the NaNs. If there are more NaNs than
% signal, then the code places NaNs for the entire signal (in that time
% bin).

%%% Inputs %%%

% LFP1, LFP2: structures with channel data
% Ex: LFP1 = load('/Volumes/SP PHD U3/AmberS/Converted Files/A7_novel1_070116/A7_novel1_070116_07-Mar-2018_CSC12.mat');
% Ex: LFP2 = load('/Volumes/SP PHD U3/AmberS/Converted Files/A7_novel1_070116/A7_novel1_070116_07-Mar-2018_CSC16.mat');
% LFP1 = load('A37_ventral hipp.mat');
% LFP2 = load('Hippocampal EEG.mat');
%
% Ex: freqRange = [2 30];
%
% Size of window of time of interest to detect coherence
% Enter 2 x 2 matrix of the beginning and end of the window of interest for
% each signal. The length of each window must be the same!
% Ex: timeWindow = [0 10; 0 10];
%
% Ex: stepWindowSize = 2;
%
% Ex: plotIT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input Checks

if LFP1.sFreq ~= LFP2.sFreq
    error('The sampling rates of your signals are different.')
end

timeStamps = LFP1.timestamps - LFP1.timestamps(1);
timeIdx1 = timeStamps >= timeWindow(1,1) & timeStamps < timeWindow(1,2);
timeStamps = LFP2.timestamps - LFP2.timestamps(1);
timeIdx2 = timeStamps >= timeWindow(2,1) & timeStamps < timeWindow(2,2);

if sum(timeIdx1) ~= sum(timeIdx2)
    error('Your singals are not the same length; please adjust your time window matrix.')
end

if timeWindow(1,2) - timeWindow(1,1) < stepWindowSize
    error('Your time window of interest is < 2 seconds making it difficult to get an accurate low frequency estimates.')
end

%% Calculation

x = LFP1.values(timeIdx1);
y = LFP2.values(timeIdx2);
newTimeStamps = 0:1/LFP1.sFreq:(timeWindow(1,2)-timeWindow(1,1))-(1/LFP1.sFreq);

stepPoints = stepWindowSize * LFP1.sFreq;
hamWin = stepPoints/4;
overlap = hamWin/2;
nFFT = [];

binSteps = (0:1000:(length(x) - stepPoints))/1000; % For coherence with smaller steps, replace 1000 (1s) with 250 (quarter of a second).

for iBin = 1:length(binSteps)
    binIdx = newTimeStamps >= binSteps(iBin) & newTimeStamps < binSteps(iBin)+stepWindowSize;
    [coherence(:,iBin),freqs] = mscohereNaN(x(binIdx),y(binIdx),hamWin,overlap,nFFT,LFP1.sFreq);
end

%% Plot

if plotIT
    
    freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);
    
    figure
    colormap(jet)
    ax(1) = subplot(211);
    plot(timeStamps(timeIdx1),[x' y'])
    ax(2) = subplot(212);
    imagesc(binSteps,1:sum(freqIdx),coherence(freqIdx,:))
    axis xy
    set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx))%,'XTick',1:5:length(binSteps),'XTickLabel',binSteps(1:5:end))
    linkaxes(ax,'x')
    
end

end