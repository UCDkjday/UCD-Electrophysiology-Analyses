function [spikeTotalAll,spikeRateAll] = countInterictalSpikes_AS(animals,mainDir,task,channels,saveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [spikeTotalAll,spikeRateAll] = countInterictalSpikes_AS(animals,mainDir,task,channels,saveDir)
%
% Counts the number of interictal spikes but also identifies the spike
% window, so it can be removed for any future spectral analyses
%
% Note: It is important to distinguish between spikes and artifact. With
% this code, we are specifically looking for spiking caused by epilepsey.
%
% INPUTS
%   animals - cell array of animal names
%   mainDir - path of where all the data are located
%   task - name of recording session you would like to count (e.g. 'novel1')
%   channels - channel in region of interest (CSC#.mat)
%   saveDir - directory where you want to save your data
%
% OUTPUTS
%   spikeTotalAll - spike counts for all animals
%   spikeRateAll - spike rate (spikes/min) for all animals
%
% Calculations
%
% 1) First, all peaks in the signal are identified using the MATLAB
% findpeaks function. With this, you can set a minimum threshold for a
% peak. This threshold is identified by the signal itself.
%     1 - Set an absolute minimum. Ex: 0.5 mV
%     2 - Calculate the std of the signal for the single channel of
%     interest -> thresholdSingleChan
%     3 - The larger of the two (absolute minimum or thresholdSingleChan)
%     is then used as the peak threshold.
%
% 2) Identifies suprathreshold peaks
%
% 3) Removes any peaks identified as artifact
%
% 4) Plots the spikes identified for each channel

% Color Guide for Plots
% black - original raw data
% blue - data identified as artifact
% green - posivite polarity spikes
% cyan - negative polarity spikes
% red - suprathreshold spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% Criteria for spike detection
numSTDsignal = 4.5; %4.5
absoluteMin = 0.5; %0.5

% Additional criteria to identify spikes
% You can set the countSpikesByPower variable to 0 if you only want to use
% a standard deviation of the signal to threshold spikes, or you can set it
% to 1 if you want to threshold by abrupt increases in power of the signal.
countSpikesByPower = 1;
numQuantile = 0.9;

% Additional criteria to identify spikes
% You can set the countSpikesBySlowWave variable to 0 if you only want to use
% a standard deviation of the signal to threshold spikes, or you can set it
% to 1 if you want to threshold by slow waves identified in the signal.
countSpikesBySlowWave = 0;
numSTDslow = 2;

% Spike Window
beforeWin = 100;
afterWin = 250;

%%  Loop through animals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create save directory if it does not exist
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

for iAnimal = 1:length(animals)
    
    tic
    
    if isempty(channels{iAnimal})
        continue
    end
    
    % Load channel to calculate a standard deviation for the entire recording
    chanFile = rdir([mainDir animals{iAnimal} '_' task '*/*' channels{iAnimal} '.mat']);
    if isempty(chanFile)
        error('The file was not found; please check your inputs and paths.')
    else
        disp(chanFile.name)
    end
    
    data = load(chanFile.name);
    sFreq = data.sFreq;
    badIntervals = data.bad_intervals;
    
    % Raw signal that has been filtered for 60 Hz line noise
    d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sFreq);
    values = filtfilt(d_60,data.values);
    %fvtool(d_60) %Look at filter reponse
    
    % Special code to handle files that are longer than they should be
    if ismember(animals{iAnimal},{'A4';'A6';'A8'})
        values = values(1:300000);
    end
    
    % Get spectral decomposition for signal to identify spikes by sharp
    % changes in power of the signal
    [~,power,~] = waveletDecomp_AS(logspace(log10(2),log10(54),30),values,sFreq,6);
    
    % Identifies artifacts from Artifact Rejection
    % NOTE: It is important to remove any peaks identified as artifact
    % because they are not spikes!
    artifactReject = zeros(1,size(values,2));

    for iArt = 1:size(badIntervals,1)
        if badIntervals(iArt,1) == 1 && badIntervals(iArt,2) == 2
            continue
        end
        if sum(sum(badIntervals > length(values)))
           badIntervals(badIntervals > length(values)) = length(values); 
        end
        artifactReject(badIntervals(iArt,1):badIntervals(iArt,2)) = 1;
    end
    artifactReject = logical(artifactReject);
    
    % Replace artifact data with Nans so artifact is not included in any calculations
    valuesNaN = values;
    valuesNaN(artifactReject) = NaN;
    powerNaN = power; 
    powerNaN(:,artifactReject) = NaN;
    
    % Determine threshold based on standard deviation
    % NOTE: Does not include artifact singal in the standard deviation
    threshSingleChan = nanstd(valuesNaN)*numSTDsignal;
    spikeThresh = max([absoluteMin threshSingleChan]);
    
    % Get the count for positive peaks
    %[posPeaks,posLocs,~,~] = findpeaks(values(iFile,:),sFreq,'MinPeakHeight',std(values)*numSTDsignal,'MinPeakDistance',1);%,'MaxPeakWidth',1);
    [posPeaks,posLocs,~,~] = findpeaks(values,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);%,'MaxPeakWidth',100);
    % Get the count for negative peaks
    [negPeaks,negLocs,~,~] = findpeaks(values*-1,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);%,'MaxPeakWidth',100);
    
    % Put location in samples instead of sec
    posLocs = round(posLocs*sFreq); negLocs = round(negLocs*sFreq);
    
    %% Identify suprathreshold peaks
    supraPeaksPos = posLocs(posPeaks > 2.8); supraPeaksNeg = negLocs(negPeaks > 2.8);
    supraSpikeVec = zeros(1,size(values,2));
    supraSpikeVec([supraPeaksPos supraPeaksNeg]) = 1;
    
    %% Count based on power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if countSpikesByPower
        
        % I chose quantile rather than std because the distribution is not
        % normal!
        normPower = powerNaN';
        normPower = bsxfun(@rdivide,normPower,max(normPower));
        aveNormPower = mean(normPower,2);
        spikeThresh = quantile(aveNormPower,numQuantile);
        %disp([num2str(iFile) ': ' num2str(spikeThresh)])
        
        % Find peaks in the average power spectrum
        [~,powLocs,~,~] = findpeaks(double(aveNormPower),sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);
        powLocs = round(powLocs*sFreq);
        
        % Match the power peaks to those peaks found in the data by threshold in the above section.
        % We do this because there is a slight mismatch in time between the two.
        combLocs = sort([posLocs negLocs]);
        if ~isempty(combLocs)
            for iClosest = 1:length(powLocs)
                [minDist,closestIdx] = min(abs(combLocs - powLocs(iClosest)));
                if minDist <= 50
                    powLocs(iClosest) = combLocs(closestIdx);
                end
            end
        end
        
        % Keep only the peaks identified by the power and signal thresholds
        [~,extraPos] = setdiff(posLocs,powLocs); [~,extraNeg] = setdiff(negLocs,powLocs);
        posPeaks(extraPos) = []; posLocs(extraPos) = [];
        negPeaks(extraNeg) = []; negLocs(extraNeg) = [];
        
    end
    
    %% Count based on the slow wave complex %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optional alternate method that is very similar to signal STD criteria.
    if countSpikesBySlowWave
        
        dLowPass = designfilt('lowpassfir','PassbandFrequency',5,'StopbandFrequency',10,'PassbandRipple',1,'StopbandAttenuation',60,'SampleRate',sFreq);
        lowPassData = filtfilt(dLowPass,values);
        %bpFilt = designfilt('bandpassfir','FilterOrder',20,'CutoffFrequency1',1,'CutoffFrequency2',10,'SampleRate',1000);
        %filtData = filtfilt(bpFilt,values(iFile,:));
        
        lowPassDataNaN = lowPassData;
        lowPassDataNaN(artifactReject) = NaN;
        spikeThresh = numSTDslow * nanstd(lowPassDataNaN);
        [posPeaksLow,posLocsLow,~,~] = findpeaks(lowPassData,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);
        [negPeaksLow,negLocsLow,~,~] = findpeaks(lowPassData*-1,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);
        
        posLocsLow = round(posLocsLow*sFreq); negLocsLow = round(negLocsLow*sFreq);
        
        % Match the power peaks to those peaks found in the data by threshold
        % in the above section.
        % We do this because there is a slight mismatch in time between the two.
        lowLocs = [posLocsLow negLocsLow];
        combLocs = sort([posLocs negLocs]);
        if ~isempty(combLocs)
            for iClosest = 1:length(lowLocs)
                [~,closestIdx] = min(abs(combLocs - lowLocs(iClosest)));
                lowLocs(iClosest) = combLocs(closestIdx);
            end
        end
        
        % Keep only the peaks identified by the power and signal thresholds
        [~,extraPos] = setdiff(posLocs,lowLocs); [~,extraNeg] = setdiff(negLocs,lowLocs);
        posPeaks(extraPos) = []; posLocs(extraPos) = [];
        negPeaks(extraNeg) = []; negLocs(extraNeg) = [];
        
    end
    
    %% Remove any peaks that are already identified as artifact
    [~,artIdxPos,~] = intersect(posLocs,find(artifactReject));
    [~,artIdxNeg,~] = intersect(negLocs,find(artifactReject));
    posPeaks(artIdxPos) = []; posLocs(artIdxPos) = [];
    negPeaks(artIdxNeg) = []; negLocs(artIdxNeg) = [];
    
    allLocs{1,1} = posLocs; allLocs{1,2} = negLocs;
    
    %% Identify spike window to be removed for future spectral analyses
    
    spikeIntervalsPos = []; spikeIntervalsNeg = [];
    
    if numel([posLocs negLocs]) ~=0
        
        artifactSpike = zeros(1,size(values,2));
        
        for i = 1:length(posPeaks)
            window = (posLocs(i) - beforeWin):(posLocs(i) + afterWin);
            window = window(window > 0 & window < size(values,2)); % Check to make sure window is within file size
            artifactSpike(window) = 1;
            spikeIntervalsPos(i,1) = window(1); spikeIntervalsPos(i,2) = window(end);
        end
        
        for i = 1:length(negPeaks)
            window = (negLocs(i) - beforeWin):(negLocs(i) + afterWin);
            window = window(window > 0 & window < size(values,2));
            artifactSpike(window) = 1;
            spikeIntervalsNeg(i,1) = window(1); spikeIntervalsNeg(i,2) = window(end);
        end
    end
    
    % Save spike information to original converted file
    spike_intervals = [spikeIntervalsPos; spikeIntervalsNeg];
    save(chanFile.name,'spike_intervals','-append')
    
    clearvars -except mainDir animals task channels saveDir num* absoluteMin countSpikesBy* *Win iAnimal values power sFreq artifactReject allLocs spikeInterval* supraSpikeVec spikeTotalAll spikeRateAll powerNaN

    %% Plot
    
    supraInChan = zeros(1,size(values,1));
    
    for iPlot = 1:size(values,1)
        
        % Plot the trace
        hFig = figure('Position',[100 100 1500 500]);
        plot(values(iPlot,:),'k') % plot 1st to set axes correctly
        ylim([min(values(:))-0.25 max(values(:))+0.25])
        hold on
        
        % Plot vertical lines indicating the location of the peak
        if ~isempty(allLocs{iPlot,1})
            vline(allLocs{iPlot,1},'--g')
            plot(allLocs{iPlot,1},max(values(:)),'v','MarkerFaceColor','g','MarkerEdgeColor','g')
        end
        
        if ~isempty(allLocs{iPlot,2})
            vline(allLocs{iPlot,2},'--c')
            plot(allLocs{iPlot,2},min(values(:)),'^','MarkerFaceColor','c','MarkerEdgeColor','c')
        end
        
        % Plot markers indicating those peaks identified as the suprathreshold peaks
        % Match the suprathreshold peaks to those peaks found in the data
        supraIdx = find(sum(supraSpikeVec));
        removeIdx = [NaN diff(supraIdx)] < 50; % removes peaks that are quite near each other
        supraIdx(removeIdx) = [];
        
        combLocs = sort([allLocs{iPlot,1} allLocs{iPlot,2}]);
        ctr = 0; % counts how many peaks were supra in that or another channel
        if ~isempty(supraIdx)
            for iClosest = 1:length(supraIdx)
                [minDist,closestIdx] = min(abs(combLocs - supraIdx(iClosest)));
                if minDist <= 50
                    foundLoc = combLocs(closestIdx);
                    if ismember(foundLoc,allLocs{iPlot,1})
                        plot(foundLoc,max(values(:)),'v','MarkerFaceColor','r','MarkerEdgeColor','r')
                        ctr = ctr + 1;
                    elseif ismember(foundLoc,allLocs{iPlot,2})
                        plot(foundLoc,min(values(:)),'^','MarkerFaceColor','r','MarkerEdgeColor','r')
                        ctr = ctr + 1;
                    end
                end
            end
        end
        %supraInChan(iPlot) = ctr;
        supraInChan(iPlot) = numel(supraIdx);
        
        % Plot again so its on top of all the other traces
        plot(values(iPlot,:),'k')
        
        % Plots spike intervals
        spikeIntervals = [spikeIntervalsNeg; spikeIntervalsPos];
        for iSpike = 1:size(spikeIntervals,1)
           plot(spikeIntervals(iSpike,1):spikeIntervals(iSpike,2),values(spikeIntervals(iSpike,1):spikeIntervals(iSpike,2)),'r')
           hold on
        end
        
        % Plot sections identified as artifact and will be removed from
        % analyses
        plot(find(artifactReject),values(iPlot,artifactReject),'.b')
        title(animals{iAnimal})
        % Plot average noramlized power
        %normPower = powerNaN(:,:,iPlot)';
        %normPower = bsxfun(@rdivide,normPower,max(normPower));
        %aveNormPower = mean(normPower,2);       
        %plot(aveNormPower,'g')
        
        % Plots low freqs
        %dLowPass = designfilt('lowpassfir','PassbandFrequency',1,'StopbandFrequency',10,'PassbandRipple',1,'StopbandAttenuation',60,'SampleRate',sFreq);
        %plot(filtfilt(dLowPass,values(iPlot,:)),'g')
        
    end
    
    %save2pdf([saveDir 'spikeDetected_' animals{iAnimal} '_' task '_' channels{iAnimal}])
    saveas(hFig,[saveDir 'spikeDetected_' animals{iAnimal} '_' task '_' channels{iAnimal} '.png'])
    waitfor(hFig)
    
    spikeTotal = sum(cellfun(@length,allLocs));
    spikeRate = spikeTotal/(sum(~artifactReject)/60000);
    
    spikeTotalAll(iAnimal) = spikeTotal;
    spikeRateAll(iAnimal) = spikeRate;
    
    %% Quantify spikes across channels
    % spikeOverlap = spikesAllChans;
    % for iRow = 1:size(allLocs,1)
    %     idx = [allLocs{iRow,1} allLocs{iRow,2}];
    %     for iCol = 1:length(idx)
    %         otherLocs = allLocs;
    %         otherLocs(iRow,:) = [];
    %         foo = cellfun(@(x) find(x < idx(iCol) + 0.05 & x > idx(iCol) - 0.05),otherLocs,'UniformOutput',false);
    %         overlap = sum(sum(cellfun(@length,foo)));
    %         spikeOverlap(iRow,round(idx(iCol)*sFreq)) = overlap + 1;
    %     end
    % end
    
    %% Save spike count for that file
    save([saveDir 'spikeTotals_' animals{iAnimal} '_' task '_' channels{iAnimal}],'values','artifactReject','spikeIntervalsPos','spikeIntervalsNeg','spikeTotal','spikeRate','supraInChan')
    
    clearvars -except mainDir animals task channels saveDir num* absoluteMin countSpikesBy* *Win spikeTotalAll spikeRateAll
    
    toc
    
end

end
