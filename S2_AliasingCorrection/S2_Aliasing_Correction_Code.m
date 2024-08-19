
% Description: Aliasing correction simulation for undersampled bSSFP profiles

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland

% 1) Profile generation based on paper:
% Plähn, N. M. J.; Poli, S.; Peper, E. S.; Açikgöz, B. C.; Kreis, R.; Ganter, C.; Bastiaansen, J. A. M. Getting the Phase Consistent:
% The Importance of Phase Description in Balanced Steady‐state Free Precession MRI of Multi‐compartment Systems. Magnetic Resonance 
% in Med 2024, mrm.30033. https://doi.org/10.1002/mrm.30033.

% 3) Used parameters: 
% M0:       polarized magnetization of the substance, i.e. PD ( 1^H proton density)
% T1:       longitudinal relaxation time
% T2:       transversal  relaxation time
% alpha:    excitation angle of RF pulse
% phi:      linear phase increment of RF excitation pulse
% TR:       repetition time of each PC-bSSFP module
% TE:       echo time of each PC-bSSFP module
% theta:    accumulated phase 
% nPC:      number of (uniformly) sampled RF phase increments 
% nvari:    variations of accumulated phase "theta" between 0 and 4*pi
% Lambda:   T1/T2 ratio

nT2      = 100;
T2t      = linspace(0.001,2,nT2);
nLambda  = 100;
Lambda   = linspace(1,60,nLambda); % starts at 1 because T1>T2 in MRI
nPCt     = 3:100; % sampling of: the number of (uniformly) sampled RF phase increments

Amount_nPC = numel(nPCt);

percent = 0.5*10^(-13);
nvari = 100;
dxi = zeros(2,Amount_nPC);
dA  = zeros(2,Amount_nPC);
dB  = zeros(2,Amount_nPC);
dz  = zeros(2,Amount_nPC);

for indnPC = 1:Amount_nPC
    nPC = nPCt(indnPC);

    for indT2 = 1:nT2
        for indLambda = 1:nLambda
            T2 = T2t(indT2);
            T1 = Lambda(indLambda)*T2;
    
            nvari = round(2*pi/nPC/deg2rad(0.5));
            if nvari<1
                nvari = 1;
            end

            nPC_infty =1000;

            for indvari = 1:nvari
                [A_inf,B_inf,z_inf,A_fin,B_fin,z_fin,Ab,Bb,zb,xib,xi_fin] = S2_Fixed_Point_validation(T1,T2,nPC,percent,indvari,nvari,nPC_infty);                    
                dxi_DFT  = xi_fin;
                dxi_bSSFP= xib;
                dA_DFT   = abs(A_inf-A_fin)/abs(A_inf)*100;
                dB_DFT   = abs(B_inf-B_fin)/abs(B_inf)*100;
                dz_DFT   = abs(z_inf-z_fin)/abs(z_inf)*100;
                dA_bSSFP = abs(A_inf-Ab)/abs(A_inf)*100;
                dB_bSSFP = abs(B_inf-Bb)/abs(B_inf)*100;
                dz_bSSFP = abs(z_inf-zb)/abs(z_inf)*100;
                % save only the MAXIMAL ERRORS to determine maximal error
                % bound -> everything is equal/below this bound
                if dxi_bSSFP>dxi(2,indnPC)
                    dxi(2,indnPC) = dxi_bSSFP;
                end
                if dA_bSSFP>dA(2,indnPC)
                    dA(2,indnPC) = dA_bSSFP;
                end
                if dB_bSSFP>dB(2,indnPC)
                    dB(2,indnPC) = dB_bSSFP;
                end
                if dz_bSSFP>dz(2,indnPC)
                    dz(2,indnPC) = dz_bSSFP;
                end
    
                if dxi_DFT>dxi(1,indnPC)
                    dxi(1,indnPC) = dxi_DFT;
                end
                if dA_DFT>dA(1,indnPC)
                    dA(1,indnPC) = dA_DFT;
                end
                if dB_DFT>dB(1,indnPC)
                    dB(1,indnPC) = dB_DFT;
                end
                if dz_DFT>dz(1,indnPC)
                    dz(1,indnPC) = dz_DFT;
                end
            end
            
        end
      
    end
    
end


%% Plot results
figure(833)
subplot(2,2,1)
loglog(nPCt,dA(2,:))
title('\Delta A in %')
subplot(2,2,2)
loglog(nPCt,dB(2,:))
title('\Delta B in %')
subplot(2,2,3)
loglog(nPCt,dz(2,:))
title('\Delta z in %')
subplot(2,2,4)
loglog(nPCt,dxi(2,:))
title('\Delta xi in %')





