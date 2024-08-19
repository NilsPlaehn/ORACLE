
function [T1_est,T2_est, theta,M0] = ORACLE_fct(bSSFPSignal,TR,alpha)
    NPhaseCycle = numel(bSSFPSignal);
    phi_tmp     = linspace(0,2*pi,NPhaseCycle+1); 
    phi         = phi_tmp(1:NPhaseCycle);      
   
    % VZ depends on the handedness of the coordinate system
    VZ = 1; 

    cm1 = conj(NPointFT(bSSFPSignal,-1.*VZ,phi));
    c0  = NPointFT(bSSFPSignal,0.*VZ,phi);
    c1  = NPointFT(bSSFPSignal,1.*VZ,phi);

    z       = c1/c0;
    theta   = angle(z);

    A        = c0;
    
    B        = cm1/z;
 
    q  = abs(A)./abs(B);

    r   = abs(z);
    
    M0 = abs(A)+abs(B);

    E2      = (1+q).*r./(r^2+q);
    T2_est = -TR./log(E2);

    Zahler  = E2-2*r+E2*r^2+E2*(1-2*E2*r+r^2)*cos(alpha);   % Definition 1
    Nenner  = E2*(1-2*E2*r+r^2)+(E2-2*r+E2*r^2)*cos(alpha); % Definition 2
    E1      = Zahler/Nenner; 
    T1_est  = -TR./log(E1);

    M0 = M0/(2*tan(alpha/2))*exp(TR/T2_est/2);
end



