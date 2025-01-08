function rfmt=RF_MT(T2c,w1,dw,lineshape)
% calculates the MT-lineshape
% last change: 2014/07/03 by PS
%w1 = gamma_2pi*x0(1,2); B1 field power in radian
%dw offset in radian Jiadi xu
if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
    step=1/1000;
    sup_theta = 0:step:pi/2;
    
    if numel(w1) == 1
        cutoff=1;
        for j=1:length(dw)
            if abs(dw(j)) >= cutoff     % the superlorentz has a pole: this avoids infinities
                                        % see Morrsion and Henkelman 1995.  units = s.  Seems weird to me. ..Daniel Gochberg  
                du=.0001;
                u=0:du:1;
                integrand2=sqrt(2/pi)*T2c./abs(3*u.^2-1) .* exp(-2*(dw(j)*T2c./abs(3*u.^2-1)).^2);
                G2(j)=w1.^2.*pi.*sum(integrand2)*du;
            else
                X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
                Y = superLorentzian_extrap_DG(T2c,w1,X,cutoff);
                G2(j)=interp1(X,Y,dw(j),'spline');
            end            
            rfmt=G2;   
        end
    else
        dwc=dw;
        Y =w1.^2.*pi.* sin(sup_theta).*sqrt(2./pi).*T2c./abs(3.*cos(sup_theta).^2-1).*exp(-2.*((dwc.*T2c)./abs(3.*cos(sup_theta).^2-1)).^2);
        rfmt=Y;
        numel(rfmt)
    end
   
elseif strcmp(lineshape,'Gaussian') %Gaussian
    rfmt=w1.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); 
elseif strcmp(lineshape,'Lorentzian') %Lorentzian
    rfmt=w1.^2*T2c./(1+(dw.*T2c).^2);
else
    error('Unknown MT-lineshape - choose SuperLorentzian, Lorentzian or Gaussian');
end;
    