function varargout = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag)
% function varargout = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag) 
% 
% Examples:
% wcut = densityCompensationForUniform3DRadial(Np, Nspoke);
% [wcut w] = densityCompensationForUniform3DRadial(Np, Nspoke);
% [wcut w] = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag);
%
% Compute density compensation weights for 3D radial acquisition. Returns
% density compensation function with correction for Nyquist limit (i.e.,
% flattened dcf outside of Nyquist limit)
% 
% INPUT:    - Np                Number of points per spoke
%           - Nspoke            Number of spokes per 3D radial acquisition
%           - oversamplingFlag  Flag to consider factor 2 oversampling
%                               (Siemens raw data)
%
% OUTPUT:   - varargout Output cell containing either wcut or {wcut w}
%               - w         Density compensation function with Nyquist condition      
%               - wcut      Density compensation function flattened out outside
%                           of Nyquist limit
%
% For computation details, see John Pauly paper on Non-Cartesian Reconstruction
%
% (c) Jerome Yerly 2013
        
  % Check input and output arguments
    narginchk(2,3); % Two input arguments
    nargoutchk(0,2); % 1 or 2 output agruments

    if nargin == 2
      % Consider no oversampling by default
        oversamplingFlag = 0;
    end
    
  % Assume a FOV
    FOV = 1;
    
  % Define idices of samples
    n = (-Np/2:1:Np/2-1).'; % Center sample is at indice N/2+1
    
  % Assume this radial sample spacing 
    dk = 1 / FOV;
    
  % Calculate density compensation by assuming dk to be 1
    w = 2*pi / (3* Nspoke) * dk^3 * (3*n.^2+1/4); 

  % Calculate area of center circle
    w(n==0) = pi / (6*Nspoke) * dk^3;
    
  % Find Nyquist limit
    nyqDiameter = sqrt( 2 * Nspoke / pi ); % See lab notebook
    
  % Correct for readout oversampling
    if oversamplingFlag
        nyqDiameter = 2 * nyqDiameter; % The factor of two considers readout oversampling
    end
    
  % Correct for oversampling case
    if (nyqDiameter > Np)
       nyqDiameter = Np; 
    end
    
  % Copy density compensation
    wcut = w;
    
  % Compute the max value before cutting (flattening out) the edges of k-space
    mask = (abs(n) <= (nyqDiameter/2));
    nyqDensity = max(abs(wcut(mask==1)));
    
  % Clip (i.e., flatten out) dcf outside of Nyquist limit
    wcut(mask==0) = nyqDensity;
    
%   % Normalize dcf
%     w = w / max(abs(w(:)));
%     wcut = wcut / max(abs(wcut(:)));

  % Return requested output
    varargout{1}    = wcut;
    if nargout == 0 || nargout == 2
        varargout{2} = w;
    end

    
    
    
%{
function varargout = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag)
% function varargout = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag) 
% 
% Examples:
% wcut = densityCompensationForUniform3DRadial(Np, Nspoke);
% [wcut w] = densityCompensationForUniform3DRadial(Np, Nspoke);
% [wcut w] = densityCompensationForUniform3DRadial(Np, Nspoke, oversamplingFlag);
%
% Compute density compensation weights for 3D radial acquisition. Returns
% density compensation function with correction for Nyquist limit (i.e.,
% flattened dcf outside of Nyquist limit)
% 
% INPUT:    - Np                Number of points per spoke
%           - Nspoke            Number of spokes per 3D radial acquisition
%           - oversamplingFlag  Flag to consider factor 2 oversampling
%                               (Siemens raw data)
%
% OUTPUT:   - varargout Output cell containing either wcut or {wcut w}
%               - w         Density compensation function with Nyquist condition      
%               - wcut      Density compensation function flattened out outside
%                           of Nyquist limit
%
% For computation details, see John Pauly paper on Non-Cartesian Reconstruction
%
% (c) Jerome Yerly 2013
        
  % Check input and output arguments
    narginchk(2,3); % Two input arguments
    nargoutchk(0,2); % 1 or 2 output agruments

    if nargin == 2
      % Consider no oversampling by default
        oversamplingFlag = 0;
    end
    
  % Assume a FOV
    FOV = 1;
    
  % Define idices of samples
    n = abs(-Np/2:1:Np/2-1).'; % Center sample is at indice N/2+1
    
  % Assume this radial sample spacing 
    dk = 1 / FOV;
    
  % Calculate density compensation by assuming dk to be 1
    w = 2*pi / (3* Nspoke) * dk^3 * (3*n.^2+1/4); 

  % Calculate area of center circle
    w(n==0) = pi / (6*Nspoke) * dk^3;
    
  % Find Nyquist limit
    nyqDiameter = sqrt( 2 * Nspoke / pi ); % See lab notebook
    
  % Correct for readout oversampling
    if oversamplingFlag
        nyqDiameter = 2 * nyqDiameter; % The factor of two considers readout oversampling
    end
    
  % Correct for oversampling case
    if (nyqDiameter > Np)
       nyqDiameter = Np; 
    end
    
  % Copy density compensation
    wcut = w;
    
  % Compute the max value before cutting (flattening out) the edges of k-space
    mask = (n <= (nyqDiameter/2));
    nyqDensity = max(abs(wcut(mask==1)));
    
  % Clip (i.e., flatten out) dcf outside of Nyquist limit
    wcut(mask==0) = nyqDensity;
    
  % Normalize dcf
    w = w / max(abs(w(:)));
    wcut = wcut / max(abs(wcut(:)));

  % Return requested output
    varargout{1}    = wcut;
    if nargout == 0 || nargout == 2
        varargout{2} = w;
    end

    
%}
    
    