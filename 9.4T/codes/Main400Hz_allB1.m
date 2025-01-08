% INSTRUCTIONS:
%       1. Load the T1map and mask from the results of Bruker_T1_RARE.m
%       2. Fit the PLOF and provide the map

clear all;
clc; 
close all;
addpath("/Your-path-to/9.4T/codes/toolbox");
datapath_pwr = 'Your-path-to/9.4T/data';
%   Load the frequency list: 83 offsets
load('path-to-data/9.4T/data/crlist.mat');

%   Load the mask of mouse brain
load('path-to-data/9.4T/data/mask.mat');

%   Load the CEST images of a B1 of 1.6uT
load('path-to-data/9.4T/data/cestimgs.mat');

%   Load the T1 map
load('path-to-data/9.4T/data/T1map.mat');

FitParam.satpwr = 1.6; % saturation power (uT)
FitParam.tsat = 2; % saturation length (second) 
FitParam.Magfield = 42.58*9.4; % 9.4 T
FitParam.CalSNR = 1;
FitParam.ifshowimage = 0;
FitParam.PeakRange = [1, 5];

%   use PLOF to calculate the CEST mappings
[FitResult] = amide_process_400Hz(ROI, imgs, fullppm, T1map, FitParam);

%  display Zamide map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
imagesc(FitResult.ZamideMap, [0, max(max(FitResult.ZamideMap))]); % 0, 4
title('Zamide map (%)')
colormap(inferno)
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'Zamide_map');
%  display Zamine map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
imagesc(FitResult.ZamineMap, [0, max(max(FitResult.ZamineMap))]); % 0, 4
title('Zamine map (%)')
colormap(inferno)
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'Zamine_map');
%   display Zguan map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
imagesc(FitResult.ZguanMap, [0, max(max(FitResult.ZguanMap))]); % 0, 4
title('Zguan map (%)')
colormap(inferno)
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'Zguan_map');
%   display M0 map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
M0 = ROI .* imgs(:, :, 2);
M0 = M0./max(max(M0));
imshow(M0, [0, 1]); % 0, 4
title('M0 map (%)')
colormap(gray(255));
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'M0_map');