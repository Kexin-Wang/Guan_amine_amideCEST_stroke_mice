function z_simulation = blochFittingBG(x, Freq, FitParam)
%BLOCHFITTINGCURVE 此处显示有关此函数的摘要
%x(1):water_T1  s
%x(2);water T2 s
%x(3); mtc T2
%x(4); mtc concentration ratio to water

% pool    = {name,            T1,      T2, lifetime,chemicalShift,concentration}
freewater = {'freewater',      x(1),     x(2),       1,        0,   110000};
cest1  = {'amide',     0.8,    0.003,     1/160,       3.5,      0};

%always set the third one as MTC
mtc  = {'mtc pool',     1    (1e-3)*x(3),     1/20,     x(5),    x(4)*110000};
pools = {freewater; cest1;mtc };

%shape pulse information
  

pulse(:,1)=FitParam.pwr*500/11.7; %1uT;
pulse(:,2)=0;
pulse(:,3)=FitParam.SatTime;
pulseCell = {pulse};
nPulseRepeat = 1;
% pulse: [strength in Hz, phase(0:x, pi/2:y, pi:-x), duration in s]
magneticField = FitParam.magneticField;       % in Tesla
gyro = FitParam.gyro;                % gyromagnetic ratio relative to 1H

Nfit=size(Freq,1);
Nreadout=(size(pools,1)*2+1);
watersignal = zeros(Nfit,2);
     
for i=1:Nfit

    frequency=Freq(i);
    y1 = 0;
    y2Mat = blochSolve_SL(pools, frequency, pulseCell, magneticField, nPulseRepeat, gyro, y1);
    % y2Mat: three dimensional matrix of magnetization vector after the pulse
    %     dim 1: [x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...]
    %     dim 2: [after pulse 1; after pulse 2; after pulse 3;...]
    %     dim 3: [after repeat 1; after repeat 2; after repeat 3;...]
    % sequence: y1 -> [pulse1, pulse2,...] -> [pulse1, pulse2,...] -> ...
    %                 \_________________________________________________/
    %                              repeat pulseRepeat times
    %     The calculation is in a rotating frame with a frequency of the pulse
    %plot result
    watersignal(i,1)=frequency;
    %for three pools. x1 x2 x3 y1 y2 y3, z1, z2 z3. so index is 7 for water
    watersignal(i,2)=y2Mat( Nreadout,end,end);

end

z_simulation = watersignal(:,2);
end

