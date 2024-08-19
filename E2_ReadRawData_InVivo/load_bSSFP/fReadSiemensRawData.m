function [twix_obj, rawData] = fReadSiemensRawData( filepathRawData )
%--------------------------------------------------------------------------
%
%   fReadRawDataSiemens     read raw data from Siemens scanner
%
%     [twix_obj rawData] = fReadSiemensRawData( filepathRawData );
%
%     INPUT:    filepathRawData Raw file path
%
%     OUTPUT:   twix_obj        Twix object containing all the header info
%               rawData         Raw data
%
%--------------------------------------------------------------------------

%   % Add ReadRawDataSiemens directory to matlab path
%     tmpdir = which('fReadSiemensRawData'); % find fReadSiemensRawData.m
% 	[tmpdir, ~] = fileparts(tmpdir);
%     addpath([tmpdir filesep 'ReadRawDataSiemens']);

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

  % READ the complete RAW DATA
  %     Raw data format = 2*Np x Nc x Ns
  %     where Np is the number of readout point, Nc is the number of
  %     channels, Ns is the total number of shots, and the factor 2 
  %     accounts for the oversampling in the readout direction.
  %
  % Order of raw data:
  %  1) Columns
  %  2) Channels/Coils
  %  3) Lines
  %  4) Partitions
  %  5) Slices
  %  6) Averages
  %  7) (Cardiac-) Phases
  %  8) Contrasts/Echoes
  %  9) Measurements
  % 10) Sets
  % 11) Segments
  % 12) Ida
  % 13) Idb
  % 14) Idc
  % 15) Idd
  % 16) Ide
    fprintf('... read raw data ...\n');
    rawData = twix_obj.image{''}; 

  % Permuting raw data to satisfy the following convention
  %     Raw data format = 2*Np x Ns x Nc
  %     where Np is the number of readout point, Nc is the number of
  %     channels, Ns is the total number of shots, and the factor 2 
  %     accounts for the oversampling in the readout direction.
  %
  % Permute data to follow convention [Nx Ny Nz Nt Nc], i.e.,
  %  1) Columns
  %  2) Lines 
  %  3) Slices
  %  4) (Cardiac-) Phases 
  %  5) Channels/Coils
  %  6) Partitions 
  %  7) Averages
  %  8) Contrasts/Echoes
  %  9) Measurements
  % 10) Sets
  % 11) Segments
  % 12) Ida
  % 13) Idb
  % 14) Idc
  % 15) Idd
  % 16) Ide
  
  % Size of the non-squeezed raw data
    dataSize = twix_obj.image.dataSize; 
    
  % "Un-squeeze" raw data
    rawData  = reshape( rawData, dataSize ); 
   
    
    fprintf('... Done!\n')
    
  % Compute elapsed time
    %toc;

