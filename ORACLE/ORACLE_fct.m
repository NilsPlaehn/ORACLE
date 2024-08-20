% Description: ORACLE framework
% This code is for research purposes only.

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland

 

function [T1_est,T2_est, theta,M0] = ORACLE_fct(bSSFPSignal,TR,alpha)
    % uniformly distributed RF phase increments
    NPhaseCycle = numel(bSSFPSignal);
    phi_tmp     = linspace(0,2*pi,NPhaseCycle+1); 
    phi         = phi_tmp(1:NPhaseCycle);      
   
    % VZ depends on the handedness of the coordinate system (Right handed coordinate system has VZ=1 and lefthanded VZ=-1)
    VZ = 1; 
    % get central modes (alternative is to pass them to this function if aliasing correction is used)
    cm1 = conj(NPointFT(bSSFPSignal,-1.*VZ,phi));
    c0  = NPointFT(bSSFPSignal,0.*VZ,phi);
    c1  = NPointFT(bSSFPSignal,1.*VZ,phi);

    % Auxilliary definitions 1
    z        = c1/c0;
    A        = c0; 
    B        = cm1/z;
    % Auxilliary definitions 2 
    x   = abs(A)./abs(B);
    r   = abs(z);   

   
    % T2 quantification
    E2      = (1+x).*r./(r^2+x);
    T2_est  = -TR./log(E2);

    % T1 quantification
    a = E2-2*r+E2*r^2;
    b = E2*(1-2*E2*r+r^2);
    E1      = (a+b*cos(alpha))/(b+a*cos(alpha)); 
    T1_est  = -TR./log(E1);
    
    % PD quantification    
    M0  = abs(A)+abs(B);
    M0  = M0/(2*tan(alpha/2))*exp(TR/T2_est/2);

    % B0 inhomogeneity/off-resonance quantification
    theta    = angle(z);
    % if df or dB desired by divide by "2*pi*TR" or "gamma*TR" respectively
    
end



