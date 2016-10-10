clear all

%% --- tracking parameters

maxMoveDist = 20; % in unit of micrometers

minIntRatio = 2.5; % minimum intensity increase in nucleus
% vs. surrounding cytoplasm

intChannel = 1; % Channel to get intensity ratios from

inputFile = 2; % 0 - new file, 1,2,{3,4} correspond to known data sets

%% --- Load raw analysis results

if inputFile == 0
    
    [sourceFile,sourceDir] = uigetfile('*.*');
    
elseif inputFile == 1
    
    sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
    sourceFile = ...
        '10_G1_Subset_fullTimeCourse.czi_AnalysisOutput.mat';

elseif inputFile == 2
    
    sourceDir = '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/';
    sourceFile = ...
        '10_G2_Subset.czi_AnalysisOutput.mat';

elseif inputFile == 3
    
    sourceDir = ...
        '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
    sourceFile = ...
        'Image_16_Subset.czi_AnalysisOutput.mat';

elseif inputFile == 4
    
    sourceDir = ...
        '/Users/hilbert/Desktop/Imaging/PCNA_SPIM_Volumes/PCNASubsets_5Feb2016/';
    sourceFile = ...
        'Image_18_Subset.czi_AnalysisOutput.mat';

end

thisPath = fullfile(sourceDir,sourceFile);

load(thisPath);

%% -- Initial plotting and evaluation

figure(1)

subplot(1,4,1)
plot(tt_vector./60,numNuc,'k-')

xlabel('Time [min]')
ylabel('N_{nuc}')


subplot(1,4,2)
plot(tt_vector./60,NN_median_vec,'k-')

xlabel('Time [min]')
ylabel('d_{NN} [\mum]')


subplot(1,4,3)
plot(tt_vector./60,cyto_vol_vec,'k-')

xlabel('Time [min]')
ylabel('V_{cyto} [\mum^3]')


ratioVec = zeros(1,numFrames);

for kk = 1:numFrames
   
   if numel(nucInt_cell{kk}{1})>0
   
       ratioVec(kk) = mean(nucInt_cell{kk}{1}./cytoInt_cell{kk}{1});
    
   else
       
      ratioVec(kk) = NaN; 
       
   end
   
end


subplot(1,4,4)
plot(tt_vector./60,ratioVec,'k-')

xlabel('Time [min]')
ylabel('Mean nuc/cyto ratio')


%% --- Nucleus tracking

trackInds = cell(1,numFrames);
trackInds{1} = 1:numNuc(1);

rvsInds = cell(1,numFrames-1);
fwdInds = cell(1,numFrames-1);

for kk = 2:numFrames
    
    if numNuc(kk)>0 && numNuc(kk-1)>0
        
        xCoords = cellfun(@(elmt)elmt(1),centroid_cell{kk-1}).';
        yCoords = cellfun(@(elmt)elmt(2),centroid_cell{kk-1}).';
        zCoords = cellfun(@(elmt)elmt(3),centroid_cell{kk-1}).';
        
        propMatrix_last = [xCoords,yCoords,zCoords];
        
        xCoords = cellfun(@(elmt)elmt(1),centroid_cell{kk}).';
        yCoords = cellfun(@(elmt)elmt(2),centroid_cell{kk}).';
        zCoords = cellfun(@(elmt)elmt(3),centroid_cell{kk}).';
        
        propMatrix_this = [xCoords,yCoords,zCoords];
        
        dist = pdist2(propMatrix_last,propMatrix_this,'seuclidean');
        
        [fwdDist,fwdInds{kk-1}] = min(dist,[],2);
        [rvsDist,rvsInds{kk-1}] = min(dist,[],1);
        
        trackInds{kk} = trackInds{kk-1}(rvsInds{kk-1});
        
    elseif numNuc(kk)>0
        
        trackInds{kk} = 1:numNuc(kk);
        
    end
    
end

% Get tracking connection length statistics, to reject jumps
moveDist = cell(1,numFrames-1);
moveDistPrctl = zeros(1,numFrames-1);
moveDistMax = zeros(1,numFrames-1);

for kk = 1:(numFrames-1)
    
    if numNuc(kk)>0 && numNuc(kk+1)>0
        
        moveDist{kk} = zeros(1,numNuc(kk+1));
        
        for nn = 1:numNuc(kk+1)
            
            currentCent = centroid_cell{kk+1}{nn};
            lastCent = centroid_cell{kk}{rvsInds{kk}(nn)};
            
            moveDist{kk}(nn) = sqrt(sum((currentCent-lastCent).^2));
            
        end
        
        moveDistPrctl(kk) = prctile(moveDist{kk},[97]);
        moveDistMax(kk) = max(moveDist{kk});
        
    end
    
end

for kk = 2:numFrames
    
    if numNuc(kk-1)>0 && numNuc(kk)>0
        
        for nn = 1:numNuc(kk)
            
            currentCent = centroid_cell{kk}{nn};
            lastCent = centroid_cell{kk-1}{rvsInds{kk-1}(nn)};
            
            % target intensity ratio check
            currentIntRatio = ...
                nucInt_cell{kk}{intChannel}(nn) ...
                ./cytoInt_cell{kk}{intChannel}(nn);
            
            
            if sqrt(sum((currentCent-lastCent).^2))>maxMoveDist ...
                    || currentIntRatio<minIntRatio
                
                rvsInds{kk-1}(nn) = NaN;
                
            end
        end
        
    end
    
end


%% --- Construct traces from tracking and trace interuptions

lastPoint = numFrames;

traces = cell(1,numNuc(lastPoint));

traceAtObjectRegister = cell(1,numFrames);

traceAtObjectRegister{numFrames} = 1:numNuc(lastPoint);

for nn = 1:numel(traces)
    
    traces{nn} = struct();
    traces{nn}.time = tt_vector(lastPoint);
    traces{nn}.timeInd = lastPoint;
    traces{nn}.origObjInd = nn;
    traces{nn}.centroid = centroid_cell{lastPoint}(nn);
    traces{nn}.volume = nuc_vol_cell{lastPoint}(nn);
    traces{nn}.intensity = ...
        cellfun(@(elmt)elmt(nn),nucInt_cell{lastPoint});
    traces{nn}.cytoIntensity = ...
        cellfun(@(elmt)elmt(nn),cytoInt_cell{lastPoint});
    
    traces{nn}.truncated = false;
    traces{nn}.joint = [];
    
    traceAtObjectRegister{numFrames}(nn) = nn;
    
end

truncatedVec = false(1,numNuc(lastPoint));
jointVec = false(1,numNuc(lastPoint));

for kk = (lastPoint-1):-1:1
    
    kk
    
    if numNuc(kk) == 0
        % Make sure all traces are truncated
        
        for tt = 1:numel(traces)
            
            if ~truncatedVec(tt) && ~jointVec(tt)
                
                traces{tt}.truncated = true;
                truncatedVec(tt) = true;
                
            end
            
        end
        
    else
        
        % positions of original objects detected in this frame, while this
        % vector holds the index of the trace associated with this original
        % object
        usedOrigInds = zeros(1,numNuc(kk));
        
        % Register that points from every original object in every frame to
        % the trace that it has been placed into
        traceAtObjectRegister{kk} = zeros(1,numNuc(kk));
        
        for tt = find(~truncatedVec & ~jointVec)
            
            if ~isempty(rvsInds{kk})
                if isfinite(rvsInds{kk}(traces{tt}.origObjInd(end)))

                    connectInd = rvsInds{kk}(traces{tt}.origObjInd(end));

                    % Check if another trace already has taken the object to
                    % connect, if yes then join this trace to that trace
                    joinInd = find(usedOrigInds==connectInd,1,'first');

                    if ~isempty(joinInd)
                        % Save time point and target trace that this trace was
                        % joint with
                        traces{tt}.joint = [kk,joinInd];
                        jointVec(tt) = true;
                    else

                        usedOrigInds(tt) = connectInd;

                        traceAtObjectRegister{kk}(connectInd) = tt;

                        traces{tt}.time = [traces{tt}.time,tt_vector(kk)];
                        traces{tt}.timeInd = [traces{tt}.timeInd,kk];
                        traces{tt}.origObjInd = ...
                            [traces{tt}.origObjInd,connectInd];
                        traces{tt}.centroid = ...
                            [traces{tt}.centroid,centroid_cell{kk}(connectInd)];
                        traces{tt}.volume = ...
                            [traces{tt}.volume,nuc_vol_cell{kk}(connectInd)];
                        traces{tt}.intensity = ...
                            [traces{tt}.intensity,...
                            cellfun(@(elmt)elmt(connectInd),nucInt_cell{kk})];
                        traces{tt}.cytoIntensity = ...
                            [traces{tt}.cytoIntensity,...
                            cellfun(@(elmt)elmt(connectInd),cytoInt_cell{kk})];
                    end
                
                else
                    
                    traces{tt}.truncated = true;
                    truncatedVec(tt) = true;
                    
                end
            end
            
        end
                        
        % Make new traces for objects whose indices were not used
        origInds = 1:numNuc(kk);
        unusedInds = setdiff(origInds,usedOrigInds);
        
        addTracesCell = cell(1,numel(unusedInds));
        truncatedVecAddition = false(1,numel(unusedInds));
        jointVecAddition = false(1,numel(unusedInds));
        
        for jj = 1:numel(unusedInds)
            
            addInd = unusedInds(jj);
            
            addTrace = struct;
            
            addTrace.time = tt_vector(kk);
            addTrace.timeInd = kk;
            addTrace.origObjInd = addInd;
            addTrace.centroid = centroid_cell{kk}(addInd);
            addTrace.volume = nuc_vol_cell{kk}(addInd);
            addTrace.intensity = ...
                cellfun(@(elmt)elmt(addInd),nucInt_cell{kk});
            addTrace.cytoIntensity = ...
                cellfun(@(elmt)elmt(addInd),cytoInt_cell{kk});
            
            addTrace.truncated = false;
            addTrace.joint = [];
            
            addTracesCell{jj} = addTrace;
                        
            traceAtObjectRegister{kk}(addInd) = numel(traces);
            
        end
        
        traces = [traces,addTracesCell];
        truncatedVec = [truncatedVec,truncatedVecAddition];
        jointVec = [jointVec,jointVecAddition];
        
    end
    
end



%% -- Overview of all traces

minTimePoints = 2; % Minimum number of time points a trace must be present
minIntRatioPeak = minIntRatio; % Minimum peak intensity ratio

figure(2)

clf

validTraceInds = ...
    cellfun(@(elmt)numel(elmt.time)>=minTimePoints,traces) ...
    & cellfun(@(elmt)max(elmt.intensity./elmt.cytoIntensity),traces) ...
    >=minIntRatioPeak;

validTraces = traces(validTraceInds);

for kk = 1:numel(validTraces)
    
    kk
    
    plot(validTraces{kk}.time./60.0,validTraces{kk}.volume,'k-')
    hold on
    
end

hold off

xlabel('Time [min]')
ylabel('V_{nuc} [\mum^3]')


%% Extract sum of nuclear volumes for different stages

% -- set parameters and context

if inputFile == 1
    
    timeWindows = ...
        {[4,14],[20,29],[32,45],[48,60],[65,77],[81,96],...
        [102,122],[124,145],[148,182]};
    
    window_labels = {'32','64','128','256','512','1K',...
        'high','oblong','sphere'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [5,6,7,8,9,10,11,12,13];
    
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

elseif inputFile == 4
    
    timeWindows = ...
        {[0,7],[12,21],[28,39],...
        [45,58],[66,86]};
    
    window_labels = {'128','256','512',...
        '1K','high'};
    
    % cell cycles corresponding to stages
    windowCellCycles = [7,8,9,10,11];
    
end


% --- execute analysis

numWindows = numel(timeWindows);

windowNumNuc = zeros(1,numWindows);
windowTotalVol = zeros(1,numWindows);
windowIndivVol = zeros(1,numWindows);
windowCytoVol = zeros(1,numWindows);

for ww = 1:numWindows
    
    ww
    
    minTime = timeWindows{ww}(1);
    maxTime = timeWindows{ww}(2);
    
    inclInd = tt_vector./60>=minTime & tt_vector./60<=maxTime;
    
    windowCytoVol(ww) = max(cyto_vol_vec(inclInd));
    
    pickTraceMask = false(size(validTraces));
    peakVols = zeros(size(validTraces));
    
    for kk = 1:numel(validTraces)
        
        thisTrace = validTraces{kk};
        
        traceTimes = thisTrace.time;
        
        traceInts = thisTrace.intensity./thisTrace.cytoIntensity;
                
        firstIntTime = traceTimes(find(traceInts>=minIntRatio,1,'last'))./60;
        lastIntTime = traceTimes(find(traceInts>=minIntRatio,1,'first'))./60;
                
        traceVols = thisTrace.volume;
        peakVols(kk) = max(traceVols);
        
        pickTraceMask(kk) =  ~(lastIntTime<minTime) && ~(firstIntTime>=maxTime);
        
    end
    
    windowNumNuc(ww) = sum(pickTraceMask);
    windowTotalVol(ww) = sum(peakVols.*pickTraceMask);
    windowIndivVol(ww) = median(peakVols(pickTraceMask));
    
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

%% --- save tracking and volume results

if inputFile == 1
    
    saveFile = 'LiveImagingVolumes_10_G1.mat';
    
elseif inputFile == 2

    saveFile = 'LiveImagingVolumes_10_G2.mat';

elseif inputFile == 3

    saveFile = 'LiveImagingVolumes_5Feb_Early.mat';

elseif inputFile == 4

    saveFile = 'LiveImagingVolumes_5Feb_Late.mat';
        
end

saveTarget = fullfile(sourceDir,saveFile);

save(saveTarget,'windowCellCycles','windowNumNuc',...
    'windowTotalVol','windowIndivVol','windowCytoVol')