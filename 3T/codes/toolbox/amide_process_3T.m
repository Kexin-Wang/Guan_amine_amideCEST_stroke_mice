function [FitResult] = amide_process_3T(mask,cestimgs,Offsets, FitParam)
%   INPUT:
%           mask: ROI of the brain
%           cestimgs: CEST images
%           Offsets: frequency list [ppm]
%           FitParam: fitting parameters
%   OUTPUT:
%           FitResult: save all the GuanCEST and amideCEST maps, etc.

%   UPDATED 2022/10/19:
%           Change the fitting peak range to (0.9, 5.0) from (1.5, 5.0) ppm
%   UPDATED 2023/02/24:
%           1. Add T1 map check 
%           2. Add the scaling of the Z-spectrum
%           3. Add CrMap
%   UPDATE 2023/03/08:
%           1. Add dZex and dZfit
%   UPDATED 2023/05/09:
%           1. Resize the cestimgs, if its size mismatched that of the mask
%   UPDATED 2023/08/09:
%           1. Add the replacement of NaN and inf in the cest images before
%           MLSVD
%   UPDATED 2023/08/11:
%           1. Add the conditional statement to avoid the case that cest
%           images are resliced by SPM and present smaller regions than the
%           masks
%   UPDATED 2023/11/10:
%           1. remove the T1 map since mouse brain is quite homogeneous
%           2. the output only contains CEST mappings of interest
%   UPDATED 2024/03/19:
%           1. Interpolate the Z-spectrum related stuff to get the results
%           in the same range as newly defined by ppm
%           2. Save it in a 64*64 array for the average of HC/ST ROI later
%   UPDATED 2024/04/11:
%           1. Get the ZamideFit and ZguanFit for each pixel

    % compute the MLSVD of the tensor
    %   Replace NaN and inf with 0 before MLSVD
    cestimgs(~isfinite(cestimgs)) = 0;
    %   All the NaN/Inf/truncated regions by SPM will be removed from the
    %   mask and the following analysis
    mask = logical(mask .* logical(squeeze(cestimgs(:,:,2))));
    [Ue, Se, ~] = mlsvd(cestimgs,[48 48 10]);
    IMGtmp_process = lmlragen(Ue, Se);


    %refine MASK image
    displayimg_process = IMGtmp_process(:,:,1);
    displayimg_process(~mask) = 0;
    NXALL = size(IMGtmp_process,1);
    NYALL = size(IMGtmp_process,2);
    for idx = 1:NXALL
       for idy = 1:NYALL
            if (mask(idx,idy))
                if (displayimg_process(idx,idy) < 0.7*mean2(displayimg_process(mask)))
                    mask(idx,idy) = false;
                end
            end
       end
    end


     CESTimg = IMGtmp_process(:,:,:)./IMGtmp_process(:,:,2);

     %remove first two M0 points
     FreqPPM = Offsets(3:end);
     CESTimg = CESTimg(:,:,3:end);
     
%      %get whole maks CEST Z-spectrum for inspection
%      for i = 1:size(FreqPPM)
%          cestinspect(:,1) = FreqPPM;
%          tmp = CESTimg(:,:,i);
%          cestinspect(i,2) = mean(tmp(mask));
%      end
     %   Frequency list
sp = -2.2; % ppm
ep = 0; % ppm
step = 0.2; % ppm
% ppm = linspace(sp, ep, n); % list of frequency in the unit of ppm
ppmn = sp: step: ep;

sp = 0.2; % ppm
ep = 4; % ppm
step = 0.1; % ppm
% ppm = linspace(sp, ep, n); % list of frequency in the unit of ppm
ppm0 = sp: step: ep;

sp = 4.2; % ppm
ep = 8.0; % ppm
step = 0.2; % ppm
% ppm = linspace(sp, ep, n); % list of frequency in the unit of ppm
ppm1 = sp: step: ep;

ppm = [ppmn,ppm0,ppm1];
ppm = flip(ppm);
ppm = ppm'; 
nrows = size(ppm,1);

     % fit each point
     idxall=0;
    
     Zamide = zeros(NXALL, NYALL);
     Zguan = zeros(NXALL, NYALL);
     Ramide = zeros(NXALL, NYALL);
     Rguan = zeros(NXALL, NYALL);
     MTbg = zeros(NXALL, NYALL);

     dZfit = zeros(NXALL, NYALL, nrows);
     dZex = zeros(NXALL, NYALL, nrows);
     Zfit = zeros(NXALL, NYALL, nrows);
     Zex = zeros(NXALL, NYALL, nrows);
     Zbg = zeros(NXALL, NYALL, nrows);
     ZamideFit = zeros(NXALL, NYALL, nrows);
     ZguanFit = zeros(NXALL, NYALL, nrows);
     
     for idx=1:NXALL
        for idy=1:NYALL
            idxall=idxall+1;

            if (mask(idx,idy))
                          
                Z_spectrum = CESTimg(idx,idy,:);
                Z_spectrum = squeeze(Z_spectrum);
       
                % B0 shift
                FitParam.WholeRange = [-2, 2];
                B0Shift = WASSR(FreqPPM, Z_spectrum, FitParam);
                
                %%   Scale the Z-spectrum to hit zero and have a standard
%                 %   amide-Z

                
                %   For debug
                if idy ==16    %   idx ==32 
                    FitParam.ifshowimage = 1; 
                    figure();
                    plot(FreqPPM, Z_spectrum);
%                     hold on;
%                     plot(xq, z_tmp);
%                     hold off;
                else
                    FitParam.ifshowimage = 0; 
                end

                    %   Scaled to [0, 1] by (Z - Zmin)./(1 - Zmin)
                    [Zmin, min_idx] = min(Z_spectrum);
                    z_tmp = (Z_spectrum - Zmin)./(1 - Zmin);
                                
                %% Fit PLOT 
                FitParam.WholeRange = [0, 8];    % CEST peak parameters
                if z_tmp(1) > 3
                    mask(idx,idy) = 0;
                    continue
                end
                [FitResult,FitParam] = PLOF(FreqPPM - B0Shift, z_tmp, FitParam);
                if size(FitResult.DeltaZspectrum, 1) > 50
                        Zamide(idx,idy) =100*FitResult.DeltaZpeak1;
                        Ramide(idx,idy) = 1000*FitResult.DeltaRpeak1;
                        Zguan(idx,idy) =100*FitResult.DeltaZpeak2;
                        Rguan(idx,idy) = 1000*FitResult.DeltaRpeak2;
                        MTbg(idx,idy) = 100 * (1 - FitResult.MTbackground);
        %                 Zcr(idx,idy) = Zguan(idx,idy) - scale.aar * Zamide(idx,idy);
                        %   (1:75) may be changed if different crlist is used
                        dZex_inter = interp1(FitResult.Offset, FitResult.DeltaZspectrum, ppm, 'spline');
                        dZex(idx,idy, :) = dZex_inter;
                        dZfit_inter = interp1(FitResult.xindex, FitResult.DeltaFitZ, ppm, 'spline');
                        dZfit(idx,idy, :) = dZfit_inter;
                        Zfit_inter = interp1(FitResult.xindex, FitResult.Curve, ppm, 'spline');
                        Zfit(idx,idy, :) = Zfit_inter;
                        Zex_inter = interp1(FitResult.Offset, FitResult.Saturation, ppm, 'spline');
                        Zex(idx,idy, :) = Zex_inter;
                        Zbg_inter = interp1(FitResult.xindex, FitResult.Background, ppm, 'spline');
                        Zbg(idx,idy, :) = Zbg_inter;
                        ZamideFit_inter = interp1(FitResult.xindex, FitResult.ZamideFit, ppm, 'spline');
                        ZamideFit(idx,idy, :) = ZamideFit_inter;
                        ZguanFit_inter = interp1(FitResult.xindex, FitResult.ZguanFit, ppm, 'spline');
                        ZguanFit(idx,idy, :) = ZguanFit_inter;
                else
                        mask(idx,idy) = false;
                end

                                               
            end
            waitbar(idxall/(NXALL*NYALL));
        end
     end
    % xindex = FitResult.xindex;
    % dZ_Offset = FitResult.Offset;
    clear FitResult;
    % FitResult.xindex = xindex;
    % FitResult.Offset = dZ_Offset;
    FitResult.ZamideMap = Zamide;
    FitResult.ZguanMap = Zguan;
    FitResult.RamideMap = Ramide;
    FitResult.RguanMap = Rguan;
    FitResult.MTbgMap = MTbg;
    FitResult.dZex = dZex;
    FitResult.dZfit = dZfit;
    FitResult.Zfit = Zfit;
    FitResult.mask = mask;
    FitResult.Zex = Zex;
    FitResult.Zbg = Zbg;
    FitResult.ZamideFit = ZamideFit;
    FitResult.ZguanFit = ZguanFit;
end