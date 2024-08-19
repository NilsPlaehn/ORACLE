function [wcut,w] = dcf3Dradial (N,NProj,filterFactor)
% ------------------------------
%%% Non Cartesian Reconstruction
% compute density compensation function for 3D radial
% adapted from Siemens code, nonCartesianFunctor.cpp
% ------------------------------
% Author: Gabriele Bonanno
% Date: November 2012
% ------------------------------


k = -N/2 : 1 : N/2-1;
k(end + 1) = N/2; 
% k = [0 : N/2-1 , -N/2 : -1]; 
w = zeros(1,length(k));


%    The factor of two considers readout oversampling
nyqDiameter = filterFactor * 2.0 * sqrt (2. * NProj / pi ); % Davide's

% nyqDiameter = filterFactor * 0.5 * sqrt ( NProj / pi ); 
% nyqDiameter = filterFactor * sqrt (2 * NProj / pi ); 

if (nyqDiameter > N)
   nyqDiameter = N; 
end

for i = 2: N    
        w(i) = ( pi * abs ( (k(i) + k(i+1))^3 - (k(i) + k(i-1))^3 ) ) / (6 * NProj);
end
% pDensCompMulMatrix[0][lEco] = 2.0 * pDensCompMulMatrix[m_iNColMeas -1][lEco] - pDensCompMulMatrix[m_iNColMeas -2][lEco];
w(1) = 2.0 * w(N) - w(N-1);

wcut = w;

% compute the max value before cutting the edges of the k-space
% maxNorm = max(w);
% maxNormcut = max(wcut);

% for (lI = m_iNColMeas/2; lI < m_iNColMeas; lI++)
for i = N/2+1 : N            
    if abs(k(i)) < (nyqDiameter / 2)
               nyqDensity = wcut(i);                
    end
end

% for i = N/2 : -1 : 1            
%     if abs(k(i)) < (nyqDiameter / 2)
%                nyqDensity2 = w(i);                
%     end
% end

% for i = N/2+1 : N
for i = 1 : N
    if abs(k(i)) > (nyqDiameter / 2)                
                wcut(i) = nyqDensity;
    end
end

% figure('name','DCF non normalized')
% plot(k(1:N),w(1:N),'-o')
% hold on, plot(k(1:N),wcut(1:N),'ro'), hold off

% add Hann filter here

% ////////////////
% // Normalization
% ////////////////

%             // Davide - find the maximum and use it for the normalization in order to deal with the normalization of the hanning filter
% 		dNormMaximum = pDensCompMulMatrix[m_iNColMeas-1];
% 		for (lI = 0; lI < m_iNColMeas; lI++)
%     		dNormMaximum = w[m_iNColMeas-1];

maxNorm = max(w);
maxNormcut = max(wcut);

w = w / maxNorm;
wcut = wcut / maxNormcut;

w = w(1:end-1);
wcut = wcut(1:end-1);
k = k(1:end-1);

% figure('name','normalized DCF')
% plot(k,w,'-o')
% hold on, plot(k,wcut,'ro'), hold off

% plot(k,repmat(nyqDensity2,[1 N+1]),'y'),hold off

%             if(m_iUseHannFilter == WIP_TRUE)
%             {
% 
%                 for (lI = 0; lI < m_iNColMeas; lI++)
%                 {
% 					
% 					dHanning = fabs (0.5 * ( 1 + cos(pi * (lI - m_iNColMeas / 2)/((double)m_iNColMeas/2)) )); //PK*4
%                     
% 					pDensCompMulMatrix[lI][lEco] *= dHanning;
% 
%                 }
%             }            

        

