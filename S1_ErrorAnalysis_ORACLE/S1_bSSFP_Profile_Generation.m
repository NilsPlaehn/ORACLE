function profile = S_bSSFP_Profile_Generation(M0,T1,T2,alpha,phi,TR,TE,theta)

% Description: generation of PC-bSSFP profiles for opposite signs

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland

% 1) Profile generation based on paper:
% Plähn, N. M. J.; Poli, S.; Peper, E. S.; Açikgöz, B. C.; Kreis, R.; Ganter, C.; Bastiaansen, J. A. M. Getting the Phase Consistent:
% The Importance of Phase Description in Balanced Steady‐state Free Precession MRI of Multi‐compartment Systems. Magnetic Resonance 
% in Med 2024, mrm.30033. https://doi.org/10.1002/mrm.30033.

% 2) Used parameters: 
% M0:      polarized magnetization of the substance, i.e. PD ( 1^H proton density)
% T1:      longitudinal relaxation time
% T2:      transversal  relaxation time
% alpha:   excitation angle of RF pulse
% phi:     linear phase increment of RF excitation pulse
% TR:      repetition time of each PC-bSSFP module
% TE:      echo time of each PC-bSSFP module
% theta:   accumulated phase 
% gamma:   gyromagnetic ratio for 1H protons
% B0:      main magnetic field strength
% dB0:     B0 inhomogeneity
% deltaCS: chemical shift


        %gamma = 2*pi*42.577*10^6; 

        E1    = exp(-TR./T1);
        E2    = exp(-TR./T2);
        %theta = -gamma*(dB0+deltaCS*B0)*TR;  
        
        a = M0.*(1-E1).*sin(alpha);
        b = 1-E1.*E2.^2+(E2.^2-E1).*cos(alpha);
        c = 2.*(E1-1).*E2.*cos(alpha./2).^2;
        profile = -1i.*a./(b+c.*cos(theta-phi)).*(1-E2.*exp(-1i.*(theta-phi))).*exp(-TE/T2).*exp(1i.*theta.*TE/TR); 
end