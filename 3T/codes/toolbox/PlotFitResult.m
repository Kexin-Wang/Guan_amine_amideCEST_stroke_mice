function PlotFitResult(FitResult,FitParam)

linecolors = lines(10);
Marker = {'-o','-v','-s','->','-*','-d','-^','-<','-p','-h','-x'};
figure();
hold on
plot(FitResult.Offset,FitResult.Saturation,'o','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(1,:));
plot(FitResult.Offset_background,FitResult.Saturation_background,'X','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(3,:));
plot(FitResult.xindex,FitResult.Curve,'-','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(2,:));
plot(FitResult.xindex,FitResult.Background,'-','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(4,:));

axis([min(FitResult.Offset),max(FitResult.Offset),min(FitResult.Saturation)*0.98,max(FitResult.Saturation)*1.02]);
% axis([1.6,3.6,0.55,0.85]);

maxSaturation = max(FitResult.Saturation);
minSaturation = min(FitResult.Saturation);
deltaSaturation = maxSaturation - minSaturation; 

% message = {
%         ['Peak Offset1 : ', num2str(FitResult.DeltaZpeak1Offset, 3)],...
%         ['Peak Offset2 : ', num2str(FitResult.DeltaZpeak2Offset, 3)],...
%         ['DeltaZpeak1 : ', num2str(FitResult.DeltaZpeak1,3)],...
%         ['DeltaZpeak2 : ', num2str(FitResult.DeltaZpeak2,3)],...
%         ['Rpeak1 : ', num2str(FitResult.Rpeak1, 3)],...
%         ['Rpeak2 : ', num2str(FitResult.Rpeak2, 3)],...
% %         ['SNR : ', num2str(FitResult.SNR,3)],...
% %         ['  '],...
%         ['MT : ', num2str(FitResult.MT,3)],...
%         ['Rsquare : ', num2str(FitResult.RsquareAll,3)]
% %         ['Noise std:',num2str(std(CurveFunction(FitResult.Coefficents, FitResult.Offset, FitParam) - FitResult.Saturation,1),3)]
%         };
 
% text(min(FitResult.Offset)+0.3*(max(FitResult.Offset)-min(FitResult.Offset)),maxSaturation-deltaSaturation*0.18,message,'FontSize',11);

FitParamMessage = evalc('FitParam');

endindex = strfind(FitParamMessage,'fields')+6;
FitParamMessage(1:endindex)=[];
text(max(FitResult.Offset),maxSaturation-deltaSaturation*0.55,FitParamMessage,'FontSize',8)

set(gca,'XDir','reverse');
box on
set (gcf,'Position',[100,100,500*1.5,400*1.5], 'color','w');
xlabel('Offset (ppm)','FontName','Times New Roman','FontSize',25,'fontweight','b');
ylabel('S/S0','FontName','Times New Roman','FontSize',25,'fontweight','b','Rotation',90);

legend('Acquired Data','Selected Background','Fitting Curve','Fitting Background','Location','SouthWest');
legend boxoff 
set(gca,'FontName','Arial','FontSize',16,'fontweight','b','LineWidth',3,'GridLineStyle','--','TickDir','in');
hold off
% 
% datapath = FitParam.datapath;
% fname = FitParam.name;
% datapath = [datapath filesep fname];
% mymkdir(datapath);
% fitsave = [datapath filesep 'FitResult.mat'];
% save(fitsave, 'FitResult');
% SaveEps(datapath, 'deltaZ')

figure();
hold on;
plot(FitResult.Offset, FitResult.DeltaZspectrum ,'o','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(1,:));
plot(FitResult.xindex,FitResult.DeltaFitZ,'-','LineWidth',3.5,'MarkerSize',11,'Color',linecolors(3,:));
xlabel('Offset (ppm)','FontName','Times New Roman','FontSize',25,'fontweight','b');
ylabel('Z','FontName','Times New Roman','FontSize',25,'fontweight','b','Rotation',90);
hold off;

end