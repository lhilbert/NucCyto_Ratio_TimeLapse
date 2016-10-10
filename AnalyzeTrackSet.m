clear all

[load_sourceFile,load_sourceDir] = uigetfile('*.*');
thisPath = fullfile(load_sourceDir,load_sourceFile);

data = load(thisPath);

rejectFlags = data.rejectFlags;
inLimitsFlags = data.inLimitsFlags;

useInds = inLimitsFlags & ~rejectFlags;

tracks = data.tracks(useInds);

numTracks = numel(tracks);

%% plotting

tt_vector = data.tt_vector;

figure(1)

clf

for tt = 1:numTracks
       
    subplot(2,1,1)

    lh = plot(tracks{tt}.time, ...
        tracks{tt}.nucIntensity{1}./tracks{tt}.cytoIntensity{1},'k-');
    set(lh,'Color',[0,0,0,0.075])
    
    hold on
    
end


set(gca,'XLim',tt_vector([1,end]))
xlabel('Time [sec]')
ylabel('N/C Pol II Ser2Phos')

for tt = 1:numTracks
        
    subplot(2,1,2)
    
    lh = plot(tracks{tt}.time, ...
        tracks{tt}.nucIntensity{2}./tracks{tt}.cytoIntensity{2},'k-',...
        'Color',[0,0,0,0.6]);
    set(lh,'Color',[0,0,0,0.075])
    
    hold on
   
end

xlabel('Time [sec]')
ylabel('N/C H3K27ac')

set(gca,'XLim',tt_vector([1,end]))