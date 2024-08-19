
function [A_inf,B_inf,z_inf,A_fin,B_fin,z_fin,Ab,Bb,zb,xib,xi_fin] = Fixed_Point_validation(T1,T2,nPC,percent,indvari,nvari,nPC_infty)
    TR      = 5/1000;  
    TE      = TR/2;
     
    theta   = 4*pi/nvari*indvari; 
    alpha   = 15.*pi/180;
    M0      = 1; 
    
    if T2>T1
        T2=T1;
    end
    
    [A_inf,B_inf,z_inf] = GetbSSFPmodes(M0,T1,T2,alpha,TR,TE,theta,nPC_infty);
    
    
    %% Get finite modes

    [c0,c1,cm1] = GetDFTmodes(M0,T1,T2,alpha,TR,TE,theta,nPC);
    A_fin = c0; 
    z_fin = c1/c0; 
    B_fin = cm1;

    xi_fin = GetXi(A_fin,B_fin,z_fin,nPC,c0,c1,cm1);

    
    
    %% Easiest iterative solving possible
    
    [Ab,Bb,zb,xib]=S2_Iterate_DTF2BSSFP_best(c0,cm1,c1,nPC,percent);
    
    

end
function xi = GetXi(A0,B0,z0,nPC,c0,c1,cm1)

    xi = sqrt((abs(A0+B0*conj(z0)^(nPC-1)-c0)^2+abs(A0*z0+B0*conj(z0)^(nPC-2)-c1)^2+...
        abs(A0*z0^(nPC-1)+B0-cm1)^2)/(abs(c0)^2+abs(cm1)^2+abs(c1)^2))*100;

end


function [c0,c1,cm1] = GetDFTmodes(M0,T1,T2,alpha,TR,TE,theta,nPC)

    phit    = linspace(0,2*pi,nPC+1);
    phi     = phit(1:nPC);
    profile = zeros(1,nPC);
    
    for indPC = 1:nPC
        profile(indPC) = S_bSSFP_profile(M0,T1,T2,alpha,phi(indPC),TR,TE,theta);
    end
    
    c0  = NPointFT(profile,0,phi);
    c1  = NPointFT(profile,1,phi);
    cm1 = NPointFT(profile,-1,phi);


end

function [A_inf,B_inf,z_inf] = GetbSSFPmodes(M0,T1,T2,alpha,TR,TE,theta,lol)
    
    nPC     = lol;
    phit    = linspace(0,2*pi,nPC+1);
    phi     = phit(1:nPC);
    
    profile = zeros(1,nPC);
    for indPC = 1:nPC
        profile(indPC) = S_bSSFP_profile(M0,T1,T2,alpha,phi(indPC),TR,TE,theta);
    end
    
    c0 = NPointFT(profile,0,phi);
    c1 = NPointFT(profile,1,phi);
    cm1 = NPointFT(profile,-1,phi);
    
    A_inf = c0; 
    z_inf = c1/c0;
    z_infc = conj(z_inf);
    B_inf = cm1;

end


function Gp = NPointFT(MatInput,order,phi)
    N = numel(MatInput);
    tmp = 0;
    for j = 1:N
        tmp = MatInput(j).*exp(1i.*order.*phi(j))+tmp;
    end   
    Gp = 1./N.*tmp;
end
function profile = S_bSSFP_profile(M0,T1,T2,alpha,phi,TR,TE,theta)
        E1    = exp(-TR./T1);
        E2    = exp(-TR./T2);
        a = M0.*(1-E1).*sin(alpha);
        b = 1-E1.*E2.^2+(E2.^2-E1).*cos(alpha);
        c = 2.*(E1-1).*E2.*cos(alpha./2).^2;
        profile = -1i.*a./(b+c.*cos(theta-phi)).*(1-E2.*exp(-1i.*(theta-phi))).*exp(-TE/T2).*exp(1i.*theta.*TE/TR); 
end





