clear all

%% --- tracking parameters

maxMoveDist = 5; % in unit of micrometers

minIntRatio = 1.15; % minimum intensity increase in nucleus
% vs. surrounding cytoplasm

minVol = 80; % minimum object volume in cubic microns
maxVol = 1800; % maximum volume

intChannel = 2; % Channel to get intensity ratios from

inputFile = 0; % 0 - new file, 1,2,{3,4} correspond to known data sets

deltat = 90; % time interval in seconds




%% --- Load raw analysis results

[sourceFile,sourceDir] = uigetfile('*.*');
    
thisPath = fullfile(sourceDir,sourceFile);

load(thisPath);

tt_vector = 0:deltat:deltat.*(numFrames-1);


%% --- Plot max projections

for ff = 1:numFrames
        
   figure(2)
   
   clf
      
   xxCoords = cellfun(@(elmt)elmt(1),centroid_cell{ff});
   yyCoords = cellfun(@(elmt)elmt(2),centroid_cell{ff});
   zzCoords = cellfun(@(elmt)elmt(3),centroid_cell{ff});
   
   subplot(2,3,1)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeY./binning).*voxelSizeY],...
       zmaxProj_cell{ff})
   
   hold on
   
   plot(xxCoords,yyCoords,'r+')
   
   xlabel('x [\mum]')
   ylabel('y [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight
   

   
   subplot(2,3,2)
    
   imagesc([0,floor(rawStackSizeY./binning).*voxelSizeY],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       xmaxProj_cell{ff}.')

   
   hold on
   
   plot(yyCoords,zzCoords,'r+')
   
   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   subplot(2,3,3)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       ymaxProj_cell{ff}.')

   
   hold on
   
   plot(xxCoords,zzCoords,'r+')
   
   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   
   
   subplot(2,3,4)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeY./binning).*voxelSizeY],...
       zmaxProj_cell{ff})
   
   xlabel('x [\mum]')
   ylabel('y [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight
   

   
   subplot(2,3,5)
    
   imagesc([0,floor(rawStackSizeY./binning).*voxelSizeY],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       (xmaxProj_cell{ff}).')

   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   
   subplot(2,3,6)
    
   imagesc([0,floor(rawStackSizeX./binning).*voxelSizeX],...
       [0,floor(rawStackSizeZ./binning).*voxelSizeZ],...
       (ymaxProj_cell{ff}).')

   xlabel('y [\mum]')
   ylabel('z [\mum]')
   
   colormap(gray)
   
   axis equal
   axis tight

   set(gca,'YDir','normal')
   
   waitforbuttonpress;
   
end