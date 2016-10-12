clear all

%% --- 

fileTargets = [1,2,3];

numTimeLapses = numel(unique(fileTargets));

cellCycles_cell = cell(1,numTimeLapses);
numNuc_cell = cell(1,numTimeLapses);
totalVol_cell = cell(1,numTimeLapses);
indivVol_cell = cell(1,numTimeLapses);
cytoVol_cell = cell(1,numTimeLapses);

for kk = 1:numTimeLapses

    kk
    
    inputFile = kk;
    
    if inputFile == 1
        % Load first timelapse results
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'TrackSet_TimeLapse_0001_G3_Subset.czi_AnalysisOutput.mat';
                
    elseif inputFile == 2
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'TrackSet_TimeLapse_0001_G2_Subset.czi_AnalysisOutput.mat';
        
    elseif inputFile == 3
        
        sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
        saveFile = 'TrackSet_TimeLapse_0001_G1_Subset.czi_AnalysisOutput.mat';
         
    end
    
    loadStruct = load(fullfile(sourceDir,saveFile));
    
    keepTrackFlags = loadStruct.inLimitsFlags & ~loadStruct.rejectFlags;
    
    includedTracks = loadStruct.tracks(keepTrackFlags);
    
    %% --- plot the time courses of number of nuclei
    
    tt_vector = loadStruct.tt_vector;
    numTimePoints = numel(tt_vector);
    
    numNuc = zeros(1,numTimePoints);
    
    parfor ll = 1:numTimePoints

%         fprintf('%3.3f percent done\n',100.*(ll-1)./(numTimePoints-1))
        
        numNuc(ll) = sum( ...
            cellfun(@(elmt)ismember(ll,elmt.timeInd),includedTracks));
        
    end
    
    
    % Extract sum of nuclear volumes for different stages
    
    % -- set parameters and context
    
    if inputFile == 1
        
        timeWindows = ...
            {[3,7.5],[17,22.5],[32,37],[45,51],[62,66],[75,79.5],...
            [90,94.5],[106,112],[126,134],[149,159],[163,173],[210,230]};
        
        window_labels = {'8','16','32','64','128','256','512','1K',...
            'hi','obl','sph','dome'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    elseif inputFile == 2
        
        timeWindows = ...
            {[3,14],[18,28],[33,43],[48,58],[64,74],[78,83],[94,98],...
            [109,113],[129,134],[155,165],[173,183],[223,233]};
        
        window_labels = {'8','16','32','64','128','256','512','1K',...
            'high','obl','sph','dome'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
    elseif inputFile == 3
        
        timeWindows = ...
            {[4,9],[18,23],[33,38],[49,52],[62,66],[79,82],[94,97],...
            [109,113],[130,133],[153,163],[175,185],[225,235]};
        
        window_labels = {'8','16','32','64','128','256','512','1K',...
            'high','obl','sph','dome'};
        
        % cell cycles corresponding to stages
        windowCellCycles = [3,4,5,6,7,8,9,10,11,12,13,14];
        
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
            
            traceTimes = thisTrace.time./60;
            
            pickTraceMask(tt) =  max(traceTimes)>=minTime ...
                && min(traceTimes)<=maxTime;
            
            if pickTraceMask(tt)
                
                thisTimeInds = thisTrace.timeInd;
                
                volStartInd = max([minInd,min(thisTimeInds)]);
                volEndInd = min([maxInd,max(thisTimeInds)]);
                
                mappedStartInd = find(thisTimeInds==volStartInd);
                mappedEndInd = find(thisTimeInds==volEndInd);
                
                volMinInd = min(mappedStartInd,mappedEndInd);
                volMaxInd = max(mappedStartInd,mappedEndInd);
                                
                traceVols = thisTrace.volume(volMinInd:volMaxInd);
                
                peakVols(tt) = max(traceVols);
                
            end
            
        end
        
        windowNumNuc(ww) = sum(pickTraceMask);
        windowTotalVol(ww) = sum(peakVols(pickTraceMask));
        windowIndivVol(ww) = mean(peakVols(pickTraceMask));
        
    end
        
    % Store results for this data set
    cellCycles_cell{kk} = windowCellCycles;
    numNuc_cell{kk} = windowNumNuc;
    totalVol_cell{kk} = windowTotalVol;
    indivVol_cell{kk} = windowIndivVol;
    cytoVol_cell{kk} = windowCytoVol;
    
end
    

%% --- Plot and save volume data

figure(1)

subplot(1,3,1)

plotStrings = {'k-o','r-s','b-^'};

for kk = 1:numTimeLapses
    
    plot(cellCycles_cell{kk},numNuc_cell{kk},...
        plotStrings{kk})
    
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