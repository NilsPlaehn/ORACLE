function [twix_obj, K] = fReadSiemensRawData_CART_MultiEcho_AM( filepathRawData )
%--------------------------------------------------------------------------
%
%     fReadRawDataSiemens     read raw data from Siemens scanner
%
%     [twix_obj rawData] = fReadSiemensRawData( filepathRawData );
%
%     INPUT:    filepathRawData Raw file path
%
%     OUTPUT:   twix_obj        Twix object containing all the header info
%               rawData         Raw data (k-space information)
%
%--------------------------------------------------------------------------
fprintf('Start reading Siemens raw data on %s\n', datestr(now));
tic;

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
rawData = twix_obj.image.unsorted(); 

% Check for multiple echoes
if twix_obj.image.NEco ~= 1
    MEFlag = true;
    fprintf('Multi-Echo (NTE = %i) acquisition detected \n', twix_obj.image.NEco);
end

% Check for Siemens line count limit
if length(twix_obj.image.Lin) >  2^16                                       
    twix_obj.image.Lin(65537:end)=1+(twix_obj.image.Lin(65537:end)-min(twix_obj.image.Lin(65537:end)));
    if length(twix_obj.image.Lin) >  2^17
        twix_obj.image.Lin(131073:end)=1+(twix_obj.image.Lin(131073:end)-min(twix_obj.image.Lin(131073:end)));
    end
end

% Sort data
if MEFlag
    K = zeros(size(rawData,1),max(twix_obj.image.Lin),max(twix_obj.image.Par),max(twix_obj.image.Eco),size(rawData,2));
    for i=1:size(rawData,3)
        K(:,twix_obj.image.Lin(i),twix_obj.image.Par(i),twix_obj.image.Eco(i),:) = rawData(:,:,i);
    end
else
    K = zeros(size(rawData,1),max(twix_obj.image.Lin),max(twix_obj.image.Par),size(rawData,2));
    for i=1:size(rawData,3)
        K(:,twix_obj.image.Lin(i),twix_obj.image.Par(i),:) = rawData(:,:,i);
    end
end

fprintf('... Done!\n')

% Compute elapsed time
toc;

end

