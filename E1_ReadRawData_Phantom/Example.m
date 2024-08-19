

% Example of how to read out raw data 

str1= 'C:\Users\Documents';
% here the RF phase increment 198° is used
str2 = '\bSSFP_raw\meas_MID00026_FID248760_PC198_FA15_TR5p00_N4_up.dat';
path = append(str1,str2);
% add the folder "load_bSSFP" to your matlab path
ref2 = Read_and_FT_raw_data_fct(path);

size(ref2)


function ref2 = Read_and_FT_raw_data_fct(path)


    filepathRawData = path; 

    
    param.RosettaTrajectoryFlag = false; 
    
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
        [twix_obj, rawData] = fReadSiemensRawData_CART(filepathRawData);
       
      % Get raw data size
        param.Np     = double(twix_obj.image.NCol);     % Number of readout point per spoke
        param.Nshot  = double(twix_obj.image.NSeg);     % Number of heartbeats
        param.Nseg   = double(twix_obj.image.NLin/param.Nshot/twix_obj.image.NEco); % Number of segments per heartbeat
        param.Nlines = double(twix_obj.image.NLin);     % Number of acquired lines
        param.Necho  = double(twix_obj.image.NEco);     % 4th dimension: Number of echoes
        param.Ncoil  = double(twix_obj.image.NCha);     % Number of coils
        param.N5     = 1;                               % 5th dimension
        param.N6     = 1;                               % 6th dimension
        param.N7     = 1;                               % 7th dimension
        param.N8     = 1;                               % 8th dimension
        param.N9     = 1;                               % 9th dimension
        param.N10    = 1;                               % 10th dimension
    
    
   
    ref2 = squeeze(ifft3c_mri(rawData)); 
end