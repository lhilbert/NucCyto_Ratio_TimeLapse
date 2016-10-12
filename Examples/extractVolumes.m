clear all

% %% --- load a reviewed track set
% 
% [load_sourceFile,load_sourceDir] = uigetfile('*.*');
% thisPath = fullfile(load_sourceDir,load_sourceFile);
% 
% data = load(thisPath);

%% --- extract volumes for specific stages

numFiles = 1;
fileTargets = [1,2,3];

numTimeLapses = numel(unique(fileTargets));

cellCycles_cell = cell(1,numTimeLapses);
numNuc_cell = cell(1,numTimeLapses);
totalVol_cell = cell(1,numTimeLapses);
indivVol_cell = cell(1,numTimeLapses);
cytoVol_cell = cell(1,numTimeLapses);

for kk = 1:numTimeLapses
    
    cellCycles_cell{kk} = [];
    numNuc_cell{kk} = [];
    totalVol_cell{kk} = [];
    indivVol_cell{kk} = {};
    cytoVol_cell{kk} = [];

end

for kk = 1:numFiles

    inputFile = kk;
    
    if inputFile == 1
        % Load first timelapse results
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'TrackSet_TimeLapse_0001_G3_Subset.czi_AnalysisOutput.mat';
                
    elseif inputFile == 2
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'LiveImagingVolumes_10_G2.mat';
        
    elseif inputFile == 3
        
        sourceDir = ...
            '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
        saveFile = 'LiveImagingVolumes_5Feb_Early.mat';
      
    end
    
    loadStruct = load(fullfile(sourceDir,saveFile));
    
    keepTrackFlags = loadStruct.inLimitsFlags & ~loadStruct.rejectFlags;
    
    includedTracks = loadStruct.tracks(keepTrackFlags);
    
    %% --- plot the time courses of number of nuclei
    
    tt_vector = loadStruct.tt_vector;
    numTimePoints = numel(tt_vector);
    
    numNuc = zeros(1,numTimePoints);
    
    parfor ll = 1:numTimePoints

        fprintf('%3.3f percent done\n',100.*(ll-1)./(numTimePoints-1))
        
        numNuc(ll) = sum( ...
            cellfun(@(elmt)ismember(ll,elmt.timeInd),includedTracks));
        
    end
    
    figure(1)
    
    subplot(1,3,1)
    
    plot(loadStruct.tt_vector./60,numNuc,'k-')
    
    xlabel('Time [min]')
    ylabel('Nuclei')
    
    % Extract sum of nuclear volumes for different stages
    
    % -- set parameters and context
    
    if inputFile == 1
        
        timeWindows = ...
            {[3,7.5],[17,22.5],[32,37],[45,51],[61,65],[75,79.5],...
            [90,94.5],[106,112],[126,134],[150,178],[199,235],[237,300]};
        
        window_labels = {'8','16','32','64','128','256','512','1K',...
            'hi','obl','sph','dome'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    elseif inputFile == 2
        
        timeWindows = ...
            {[3,14],[18,28],[33,43],[48,58],[64,76],...
            [81,100],[108,138]};
        
        window_labels = {'64','128','256','512','1K',...
            'high','oblong'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [6,7,8,9,10,11,12];
        
    elseif inputFile == 3
        % still needs adjustment
        
        timeWindows = ...
            {[3,8],[17,24],[30,39]};
        
        window_labels = {'16','32','64'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [4,5,6];
        
    end

    
    % --- execute analysis
    
    numWindows = numel(timeWindows);
    
    windowNumNuc = zeros(1,numWindows);
    windowTotalVol = zeros(1,numWindows);
    windowIndivVol = zeros(1,numWindows);
    windowCytoVol = zeros(1,numWindows);
    
    for ww = 1:numWindows
                
        minTime = timeWindows{ww}(1);
        maxTime = timeWindows{ww}(2);
        
        minInd = find(tt_vector./60>=minTime,1,'first');
        maxInd = find(tt_vector./60<=maxTime,1,'last');
        inclInd = minInd:1:maxInd;
        
        windowCytoVol(ww) = max(loadStruct.cyto_vol_vec(inclInd));
        
        pickTraceMask = false(size(includedTracks));
        peakVols = zeros(size(includedTracks));
        
        for tt = 1:numel(includedTracks)
            
            thisTrace = includedTracks{tt};
            
            traceTimes = thisTrace.time;
            
            pickTraceMask(tt) =  min(traceTimes)>=minTime ...
                && max(traceTimes)<=maxTime;

            if pickTraceMask(tt)
                
                thisTimeInds = thisTrace.timeInd;
                
                volStartInd = max([minInd,min(thisTimeInds)]);
                volEndInd = min([maxInd,max(thisTimeInds)]);
                
                traceVols = thisTrace.volume(volStartInd,volEndInd);
                peakVols(kk) = max(traceVols);
                
            end
            
        end
        
        windowNumNuc(ww) = sum(pickTraceMask);
        windowTotalVol(ww) = sum(peakVols.*pickTraceMask);
        windowIndivVol(ww) = median(peakVols(pickTraceMask));
        
    end
    
    
    subplot(1,3,2)
    
    plot(windowCellCycles,windowNumNuc,'k-')
    
    xlabel('Stage')
    ylabel('Nuclei')
    
end
    

%% --- Total cytoplasmic volume

[cytoVolPeaks,peakInds] = findpeaks(cyto_vol_vec,'MinPeakDistance',6);

figure(3)

plot(tt_vector(peakInds)./60,cytoVolPeaks,'ko')
hold on
plot(tt_vector./60,cyto_vol_vec,'b-')
hold off

ylim([0,inf])

xlabel('Time [min]')
ylabel('V_{cyto} [\mum^3]')

cytoVol = mean(cytoVolPeaks(1:5));

    
    
        cyto_vol = loadStruct.cyto_vol_vec;

        
    
%     cellCycles_cell{fileTargets(kk)} = ...
%         [cellCycles_cell{fileTargets(kk)}, ...
%         loadStruct.windowCellCycles];
%     numNuc_cell{fileTargets(kk)} = ...
%         [numNuc_cell{fileTargets(kk)}, ...
%         loadStruct.windowNumNuc];
%     totalVol_cell{fileTargets(kk)} = ...
%         [totalVol_cell{fileTargets(kk)},...
%         loadStruct.windowTotalVol];
%     indivVol_cell{fileTargets(kk)} = ...
%         {indivVol_cell{fileTargets(kk)}, ...
%         loadStruct.windowIndivVol};
%     cytoVol_cell{fileTargets(kk)} = ...
%         [cytoVol_cell{fileTargets(kk)},...
%         loadStruct.windowCytoVol];
    
end
    

figure(1)

subplot(1,3,1)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},numNuc_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('Nuclei count')




subplot(2,3,2)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},totalVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('V_{nuc} [\mum^3]')


subplot(2,3,5)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},cytoVol_cell{kk},plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('V_{cyto} [\mum^3]')






subplot(1,3,3)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    refVol = mean(cytoVol_cell{kk}...
        (cellCycles_cell{kk}>=8 & cellCycles_cell{kk}<=10));
    
    plot(cellCycles_cell{kk},...
        totalVol_cell{kk}./refVol,...
        plotStrings{kk})
    
    hold on

end

hold off

xlabel('Cell cycle')
ylabel('%V_{nuc}')