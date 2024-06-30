function OUT = power_AS(LFP,freqs,window,PLOT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function OUT = peakFreqPower_AS(LFP,freqs,window,PLOT)
%
% INPUTS
%   LFP: LFP structure or a file name containing the signal information
%   freqs: vector of frequencies of interest
%       Ex: freqs = logspace(log10(3),log10(54),24);
%   window: time window of interest you want the mean peak frequency in sec
%       Note: [] defaults to whole recording
%       Ex: window = [0 7]; -> 0 to 7 seconds
%   PLOT: plots the raw signal, power spectrum, and peak freq power
%       Ex: PLOT = 0; -> does not output a figure
%
% OUTPUT
%   OUT: structure containing many calculated values (see code below)
%
% Example: peakFreqPower_AS([],[6 10],[10 20],0)
% Calculate peak power for a GUI selected LFP for frequencies 6 and 10 hz
% at times 10 seconds to 20 seconds
%
% Example: peakFreqPower_AS([],[20 50],[],1)
% Calculate peak power for a GUI selected LFP for frequencies 20 to 50 hz
% for the whole recording and plot it!
%
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check and modify inputs
if isempty(LFP)
    [LFP,path] = uigetfile('','Select LFP');
    LFP = fullfile(path,LFP);
    LFP = load(LFP);
elseif ischar(LFP)
    LFP = load(LFP);
end

% Create the timestamps relative to start
timeStamps = LFP.timestamps - LFP.timestamps(1);
if nargin < 3 || isempty(window)
    window = [0 timeStamps(end)];
end

% Defaults to no output figure
if nargin < 4
    PLOT = 0;
end

%% Extract power of signal in LFP

% Handle bad data by interpolation (better handles edge effects from fft)
[interpolatedData,badIdx] = interpolateBadData(LFP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE BUFFER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buffer the EEG data to ensure no edges artifacts at beginning and end
[bufferedData,bufferBins] = addMirroredBuffers(interpolatedData,freqs,LFP.sFreq);

% Extract power from raw LFP data using wavelet decomposition over
% range of frequencies specifiefd by freqs
waveletNum = 6;
[phase,power,~] = waveletDecomp_AS(freqs,bufferedData,LFP.sFreq,waveletNum);

% Trim the buffers
power = power(:,(bufferBins + 1):(length(power) - bufferBins));
phase = phase(:,(bufferBins + 1):(length(phase) - bufferBins));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE BUFFER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

artifactReject = zeros(1,length(interpolatedData));
artifactReject(badIdx) = 1; artifactReject = logical(artifactReject);

%% Mean peak freq and power calculations

% Select the correct time window
tsIdx = timeStamps >= window(1) & timeStamps < window(2);

% NaN artifact time points
power(:,badIdx) = NaN;
phase(:,badIdx) = NaN;

% Calculate mean frequency in theta range
[~,thetaIdx1] = min(abs(freqs - 6));
[~,thetaIdx2] = min(abs(freqs - 10));
[peakPower,peakIdx] = max(power(thetaIdx1:thetaIdx2,tsIdx),[],'includenan'); % find the maximum of the power in each window
thetaFreqs = freqs(thetaIdx1:thetaIdx2);
peakFreq = thetaFreqs(peakIdx); % convert those indicies to actual frequencies
peakFreq(artifactReject(tsIdx)) = NaN;
peakPower(artifactReject(tsIdx)) = NaN;
peakIdx(artifactReject(tsIdx)) = NaN;

% Condense output
OUT.mean_peak_freq = nanmean(peakFreq);
OUT.median_peak_freq = nanmedian(peakFreq);
OUT.std_peak_freq = nanstd(peakFreq);
OUT.mean_peak_power = nanmean(peakPower);
OUT.std_peak_power = nanstd(peakPower);
OUT.frequencies = freqs;
OUT.power = power;
OUT.phase = phase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot

if PLOT

    % Spectrogram with peak frequency plotted in red
    figure('Position',[100 100 1500 1000])
    ax(1) = subplot(211);
    rawData = LFP.values; rawData(artifactReject) = NaN;
    plot(timeStamps(tsIdx),rawData(tsIdx))
    title('Raw Data')
    ax(2) = subplot(212);
    imagesc(timeStamps(tsIdx),1:length(freqs),power(:,tsIdx))
    %contourf(timeStamps(tsIdx),freqs,power(:,tsIdx),20,'linecolor','none')
    axis xy
    hold on
    %plot(timeStamps(tsIdx),peakIdx,'.r')
    set(gca,'YTick',1:5:length(freqs),'YTickLabel',freqs(1:5:end))
    %set(gca,'YTick',round(logspace(log10(freqs(1)),log10(freqs(end)),10)*100)/100,'yscale','log','clim',[-90 0])
    title(['Mean peak frequency is ' num2str(OUT.mean_peak_freq) 'Hz'])
    xlabel('Time (s)')
    ylabel('Freq (Hz)')
    h = colorbar;
    ylabel(h,'dB')
    linkaxes(ax,'x')
    
    % Histogram of peak frequencies over time window of interest
    %figure('Position',[100 100 1500 1000])
    %histogram(peakFreq,length(freqs))
    %set(gca,'XTick',1:length(freqs),'XTickLabel',freqs,'XTickLabelRotation',315)
    
end

end

