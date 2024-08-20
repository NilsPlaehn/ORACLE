% Description: Coil combination code for PC bSSFP profile generation

% This code is for research purposes only.

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland


% Coil combination code is not published. 
% If questions arise please contact: nils.plaehn@students.unibe.ch




clear all


%% load the data from the scanner
load('Mat_Coils.mat');

% create an array where the x,y,z voxels with Nx Ny Nz matrix size is in
% the first 3 entries
% in the 4th entry is the number of sampled RF phase increments
% in the 5th entry is the number of coils

%% Dimensions
Nx    = size(Mat,1);
Ny    = size(Mat,2);
Nz    = size(Mat,3);
nPC   = size(Mat,4);
Ncoil = size(Mat,5);


%% PC bSSFP coil combination (not published - if question please contact me)
 
bSSFP = zeros(Nx,Ny,Nz,nPC);

for indCoil=1:Ncoil 
    tmp = squeeze(Mat(:,:,:,:,indCoil));
    CS = sum(tmp,4); % complex sum (CS) = center of ellipse as a reference we rotate all ellipses to, to get a "stack of pancakes"
    for indPC = 1:nPC
        % rotate all elliptical profiles to the real axis such that the CS
        % is real
        tmp(:,:,:,indPC) = squeeze(tmp(:,:,:,indPC)).*exp(-1i.*angle(CS));
    end
    % complex sum of different elliptical profiles for different coils
    bSSFP = bSSFP+tmp;
end

%% SF
str_save = append('Pancake.mat');
save(str_save,'bSSFP');
