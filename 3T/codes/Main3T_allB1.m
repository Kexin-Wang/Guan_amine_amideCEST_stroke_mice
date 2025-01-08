% INSTRUCTIONS:
%       1. Load CEST images of one B1 value, mask, frequency list
%       2. Fit the PLOF and estimate the maps of amideCEST and GuanCEST

clear all;
clc; 
close all;

addpath("path-to-toolbox/toolbox");
datapath_pwr = 'path-to-data/3T/data'; %    For saving the plots

%   Load the frequency list: 73 offsets
load('path-to-data/crlist.mat');

%   Load the mask of mouse brain
load('path-to-data/mask.mat');

%   Load the CEST images of a B1 of 0.8uT
load('path-to-data/cestimgs.mat');


FitParam.tsat = 2; % saturation length (second) 
FitParam.Magfield = 42.58*3; % 3 T
FitParam.R1 = 1/1.205;
FitParam.CalSNR = 1;
FitParam.ifshowimage = 0;
FitParam.PeakRange = [1, 5];
FitParam.satpwr = 0.8;
imgs = squeeze(imgs);   %   64*64*73

%   use PLOF to calculate the CEST mappings
[FitResult] = amide_process_3T(ROI, imgs, fullppm, FitParam);

%   display amideCEST map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
imagesc(FitResult.ZamideMap, [0, min(6.5, max(max(FitResult.ZamideMap)))]); % 0, 4
title('Zamide map (%)')
colormap(inferno)
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'Zamide_map');

%   display GuanCEST map
figure();
set(gca,'Position',[0.05 0.05 0.9 0.9]);
imagesc(FitResult.ZguanMap, [0, min(6, max(max(FitResult.ZguanMap)))]); % 0, 4
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
M0 = FitResult.mask .* imgs(:, :, 2);
M0 = M0./max(max(M0));
imshow(M0, [0, 1]); % 0, 4
title('M0 map (%)')
colormap(gray(255));
colorbar('location','Eastoutside','FontSize', 18)
set(gca,'dataAspectRatio', [1 1 1]);
axis off
hold off
SaveEps(datapath_pwr, 'M0_map');
