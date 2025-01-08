function dZ = blochFittingCurve(x, Freq, FitParam)

%x(1):water_T1  s
%x(2);water T2 s
%x(3); amide exchange rate
%x(4); amide concentration mM
%x(5); mtc T2 ms
%x(6); mtc ratio to water

BGpara = FitParam.BGpara;
% pool    = {name,            T1,      T2,      lifetime,chemicalShift,     concentration}
freewater = {'freewater',      BGpara(1),     BGpara(2),       1,        0,   110000};
cest1  = {'guan',     2,    (1e-3)*60,     1/x(1),      2.,      20};
%always set the third one as MTC
mtc  = {'mtc pool',     1,    (1e-3)*BGpara(3),     1/20,     BGpara(5),    BGpara(4)*110000};
pools = {freewater; cest1;mtc };
pulse(:,1) = FitParam.pwr*500/11.7; %1uT;
pulse(:,2) = 0;
pulse(:,3) = FitParam.SatTime;
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

%	For the GAMT, the Guan pool is mainly arginine pool
% pool = {name,   T1,      T2,           lifetime,       chemicalShift,     concentration}
cest1  = {'arginine',     2,    (1e-3)*x(4),     1/x(3),          2.,                         x(2)};
pools = {freewater; cest1;mtc };
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

zback = watersignal(:,2);
dZ = zback - z_simulation;
end

