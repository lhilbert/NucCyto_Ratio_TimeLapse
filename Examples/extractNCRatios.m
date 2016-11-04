clear all

%% --- prepare containers

fileTargets = [1];

numTimeLapses = numel(unique(fileTargets));

tt_vec_cell = cell(1,numTimeLapses);
numNuc_vec_cell = cell(1,numTimeLapses);

cellCycles_cell = cell(1,numTimeLapses);
numNuc_cell = cell(1,numTimeLapses);

volumes_cell = cell(1,numTimeLapses);
NCratio_cell = cell(1,numTimeLapses);

timeWindow_cell = cell(1,numTimeLapses);

%% --- prepare time windows

for kk = 1:numTimeLapses
    
   if kk == 1
        
        timeWindows = ...
            {[20,30],[35,50],[53,63],[66,77],[82,92],[96,110],...
            [112,129],[130,150]};
        
        % cell cycles corresponding to stages
        windowCellCycles = [2,3,4,5,6,7,8,9];
        
    end

    timeWindow_cell{kk} = timeWindows; 
    cellCycles_cell{kk} = windowCellCycles;
    
end


%% --- extraction of volumes

for kk = 1:numTimeLapses

    kk
    
    inputFile = kk;
    
    if inputFile == 1
        % Load first timelapse results
        
        sourceDir = '/Users/hilbert/Desktop/Imaging';
        saveFile = 'TrackSet_Yuko_EmbryoFabData_AnalysisOutput.mat';
        
    end
    
    loadStruct = load(fullfile(sourceDir,saveFile));
    
    keepTrackFlags = loadStruct.inLimitsFlags & ~loadStruct.rejectFlags;
    
    includedTracks = loadStruct.tracks(keepTrackFlags);
    
    numChannels = numel(includedTracks{1}.cytoIntensity);
    
    %% --- plot the time courses of number of nuclei
    
    tt_vector = loadStruct.tt_vector;
    numTimePoints = numel(tt_vector);
    
    numNuc = zeros(1,numTimePoints);
    
    parfor ll = 1:numTimePoints
        
        numNuc(ll) = sum( ...
            cellfun(@(elmt)ismember(ll,elmt.timeInd),includedTracks));
        
    end
    
    tt_vec_cell{kk} = tt_vector;
    numNuc_vec_cell{kk} = numNuc;
    
    % Extract sum of nuclear volumes for different stages
    
    % -- set parameters and context
    
    timeWindows = timeWindow_cell{kk};
    
    % --- execute analysis
    
    numWindows = numel(timeWindows);
    
    windowNumNuc = zeros(1,numWindows);
    windowSingleNucleiVolumes = cell(1,numWindows);
    windowSingleNucleiNCRatios = cell(1,numWindows);
        
    for ww = 1:numWindows
                
        minTime = timeWindows{ww}(1);
        maxTime = timeWindows{ww}(2);
        
        minInd = find(tt_vector./60>=minTime,1,'first');
        maxInd = find(tt_vector./60<=maxTime,1,'last');
        inclInd = minInd:1:maxInd;
                
        pickTraceMask = false(1,numel(includedTracks));
        peakVols = zeros(1,numel(includedTracks));
        
        peakNCRatios = zeros(numChannels,numel(includedTracks));
        
        for tt = 1:numel(includedTracks)
            
            thisTrace = includedTracks{tt};
            
            traceTimes = thisTrace.time./60;
            
            pickTraceMask(tt) =  max(traceTimes)>=minTime ...
                && min(traceTimes)<=maxTime;
            
            if pickTraceMask(tt)
                
                thisTimeInds = thisTrace.timeInd;
                
                trackStartInd = max([minInd,min(thisTimeInds)]);
                trackEndInd = min([maxInd,max(thisTimeInds)]);
                
                mappedStartInd = find(thisTimeInds==trackStartInd);
                mappedEndInd = find(thisTimeInds==trackEndInd);
                
                actualMinInd = min(mappedStartInd,mappedEndInd);
                actualMaxInd = max(mappedStartInd,mappedEndInd);
                                
                traceVols = thisTrace.volume(actualMinInd:actualMaxInd);
                peakVols(tt) = max(traceVols);
                
                for cc = 1:numChannels
                    cytoInts = thisTrace.cytoIntensity{cc}(...
                        actualMinInd:actualMaxInd);
                    nucInts = thisTrace.nucIntensity{cc}(...
                        actualMinInd:actualMaxInd);
                    peakNCRatios(cc,tt) = max(nucInts./cytoInts);
                end
                
            end
            
        end
        
        windowNumNuc(ww) = sum(pickTraceMask);
        windowSingleNucleiVolumes{ww} = peakVols(pickTraceMask);
        windowSingleNucleiNCRatios{ww} = peakNCRatios(:,pickTraceMask);
        
    end
    
    % Store results for this data set
    numNuc_cell{kk} = windowNumNuc;    
    NCratio_cell{kk} = windowSingleNucleiNCRatios;
    volumes_cell{kk} = windowSingleNucleiVolumes;
    
end
    

%% --- Save volume data

save([sourceDir,'LiveImagingVolumes'],...
    'numTimeLapses','numNuc_cell',...
    'volumes_cell','NCratio_cell',...
    'timeWindow_cell','cellCycles_cell')

%% --- Plot volume data

plotStrings = {'k-o','r-s','b-^'};


for kk = 1:numTimeLapses
    
    numStages = numel(timeWindow_cell{kk});
    
    for ww = 1:numStages
    
        plotColor = [0.7,0.7,0.7].*(ww-1)./(numStages-1);
        
        plot(volumes_cell{kk}{ww}(:),NCratio_cell{kk}{ww}(1,:),...
            'ko','MarkerEdgeColor',plotColor)
    
        hold on
        
    end
    
end

hold off