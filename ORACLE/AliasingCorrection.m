% Description: Aliasing correction

% This code is for research purposes only.

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland


function [Ab,Bb,zb,xib]=AliasingCorrection(c0,cm1,c1,nPC,percent)

    rho = c1/c0;
    chi = cm1/c0;

    ntmp = 10^5;
    zvec   = zeros(1,ntmp);
    zvec(1)=rho;

    xi1 = zeros(1,ntmp);    
    k   = 1; 
    xib = 10000; 


    while 1==1
        
        z = Getz(rho,chi,zvec(k),nPC);
        
        zvec(k+1) = z;


        % calculate A, B from z
        r = abs(z);        
        A = (c0-cm1*conj(z)^(nPC-1))/(1-r^(2*nPC-2));
        B = (cm1-c0*z^(nPC-1))/(1-r^(2*nPC-2));
        % calculate RMSE in % for break up condition
        xitmp  = GetXi(A,B,zvec(k+1),nPC,c0,c1,cm1);
        xi1(k) = xitmp; 
        if xitmp<xib
            xib = xitmp;
            Ab  = A;
            Bb  = B;
            zb  = z;
        end
        if or(xi1(k)<percent,k>numel(zvec))
            break; 
        end

        k = k+1;
    end

    Ab = Ab*(1-zb^nPC);
    Bb = Bb*(1-conj(zb)^(nPC));
    

end

function z = Getz(rho,chi,z,n)
    z = rho-(chi-z^(n-1))/(1-abs(z)^(2*n-2))*(1-abs(z)^(2))*conj(z)^(n-2);
end


function xi = GetXi(A0,B0,z0,nPC,c0,c1,cm1)

    xi = sqrt((abs(A0+B0*conj(z0)^(nPC-1)-c0)^2+abs(A0*z0+B0*conj(z0)^(nPC-2)-c1)^2+...
        abs(A0*z0^(nPC-1)+B0-cm1)^2)/(abs(c0)^2+abs(cm1)^2+abs(c1)^2))*100;

end

