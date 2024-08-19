% Example of how to read out raw data 
clear all
str1= 'C:\Users\Documents';
% here the RF phase increment 198° is used
str2 = 'meas_MID00021_FID239867_A3_PCBSSFP_Auto_180_V02.dat';
path = append(str1,str2);
% add the folder "load_bSSFP" to your matlab path

rawData= F_BSSFP_FourierTransform(path);

size(rawData)

%% Visualize Data

indz = 10; 
indCoil = 10;
close all
for indPC = 1:21 % RF phase increment 21 is the same as 1 : check for for B0 drift, movement etc
    figure(indPC)
    imagesc(squeeze(abs(rawData(:,:,indz,indPC,indCoil))))
    axis image
    title(indPC)
    %pause(0.5)
end

function rawData= F_BSSFP_FourierTransform(path)

    
    
        filepathRawData = path; 
        %% PART 1
        %%% Load and read 3D radial data - main
        
        param.RosettaTrajectoryFlag = false; % set to phyllotaxis by default
        
        %--------------------------------------------------------------------------
        % Miscellaneous parameters
        %--------------------------------------------------------------------------
        param.dispFlag      = true;
        param.flagSelfNav   = true;
        param.dcf_cutFlag   = true;
        
        %--------------------------------------------------------------------------
        % Load data
        %--------------------------------------------------------------------------
          % Load twix object and read raw data in format [Np Nseg*Nshot Necho Ncoil]
           % [twix_obj, rawData] = fReadSiemensRawData_CART(filepathRawData);
    
    
    %% -------------------------------------------------------------------
    %% -------------------------------------------------------------------     
    %% -------------------------------------------------------------------
    %% -------------------------------------------------------------------
    %% -------------------------------------------------------------------     
    %% -------------------------------------------------------------------
    
    
           fprintf('Start reading Siemens raw data on %s\n', datestr(now));
    %tic;
    
    % Reading RAW DATA header (raw data are actually read only when needed);
    fprintf('... read raw data header ...\n');
    twix_obj_multi = mapVBVD(filepathRawData);
    if iscell(twix_obj_multi)
        twix_obj = twix_obj_multi{2};
    else
        twix_obj = twix_obj_multi;
    end
    twix_obj.image.flagIgnoreSeg = true; %(essential to read the data correctly)
    

    fprintf('... read raw data ...\n');
    rawData = twix_obj.image.unsorted(); % Cartesian
    
    
    
    A_Lin = twix_obj.image.Lin(:);
    A_Par = twix_obj.image.Par(:);
    A_Set = twix_obj.image.Set(:);
    
    nVec  = numel(A_Lin);
    
    indtmp = 1; 
    for k =1:nVec
        if A_Lin(k)>65000
            A_Lin(k) = A_Lin(k)-65536;
        end
        if k>1
            if A_Par(k)<A_Par(k-1)
                indtmp = indtmp+1;
            end
            A_Set(k) = indtmp;
        end
    end
    
    
    
    figure(1)
    plot(A_Lin)
    figure(2)
    plot(A_Par)
    figure(3)
    plot(A_Set)
    
    max(A_Lin)
    max(A_Par)
    max(A_Set)
    
    %% sorting script Adele
    % TRAJ = single( [twix_obj.image.Lin(:), twix_obj.image.Par(:), twix_obj.image.Set(:)] );
    
    TRAJ = single( [A_Lin, A_Par, A_Set] );
    tempK = zeros(size(rawData, 1), size(rawData, 2), max(TRAJ(:,1)), max(TRAJ(:,2)), max(TRAJ(:,3)));%,'single');
    for i = 1:size(rawData, 3)
        tempK(:, :, TRAJ(i,1), TRAJ(i,2), TRAJ(i,3)) = squeeze(rawData(:,:,i));
    end
    % tempK format: nCol (Nx) x nCoils x nLin (Ny) x nPar (Nz)
    % reformat to Nx x Ny x Nz x Ncoils
    tempK = permute(tempK, [1 3 4 5 2 ]);
    rawData = tempK;
    
    fprintf('... Done!\n')
    
    %% -------------------------------------------------------------------

    
    rawData = squeeze(ifft3c_mri(rawData)); % ohne coil compression
    
%     %% SAff
%     Nx             = size(rawData,1);
%     Ny             = size(rawData,2);
%     Nz             = size(rawData,3);
%     nPC            = size(rawData,4);
%     Ncoil          = size(rawData,5);
% 
%     phi1 = linspace(0,360,nPC+1);
%     phi  = phi1(1:nPC);
%     
%     
%     cd(str_Metafolder)
%     mkdir(str_SavedBSSFP);
%     cd(currentFolder);
%     for indPC = 1:nPC
%         cd(str_SavedBSSFP)
%         for indz = 1:Nz
%             strname = append('PC',num2str(phi(indPC)),'_Slice',num2str(indz),'.mat');
%             varname = squeeze(rawData(:,:,indz,indPC,:)); 
%             save(strname,'varname');
%         end
%         cd(currentFolder);
%         indPC
%     end
end



















