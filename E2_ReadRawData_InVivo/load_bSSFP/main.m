%%%%% Load and read 3D radial data - main
%%%%% ALCM 10.03.2022



% filepathRawData     = '/mnt/MountFilearc/RAD/CIBM_CVMR/nD-Repository/_project_5D_FISS_/PRISMA_20190218_3X5D_1p4mm_1xSN_1p1mm_Aurelien/meas_MID00371_FID84159_SN_LIBR_192_T2P40_a23_d1800_f680_s20_tra24_513.dat'; % full path to rawdata (.dat) here
filepathRawData     = 'C:\Users\evasp\Documents\Projects\GRASP\3D_rad_data\meas_MID00020_FID137231_3Drad_22x110.dat'

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
    [twix_obj, rawData] = fReadSiemensRawData(filepathRawData);

   
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

%--------------------------------------------------------------------------
% Compute trajectory and density compensation weights of original data
%--------------------------------------------------------------------------   
  % Choose k-space trajectory
    if ~isfield(param, 'RosettaTrajectoryFlag')
        button = questdlg('Choose k-space trajectory','K-space trajectory selection','Phyllotaxis','Rosetta','Phyllotaxis');
        param.RosettaTrajectoryFlag = strcmp(button,'Rosetta');
    else
        param.RosettaTrajectoryFlag = false;
        disp('Phyllotaxis trajectory selected by default!');
    end
    
  % Compute phyllotaxis trajectory
    [kx, ky, kz] = computePhyllotaxis (param.Np, param.Nseg, param.Nshot, param.flagSelfNav, param.dispFlag, param.RosettaTrajectoryFlag);
    Kstatic   = cat(4, kx, ky, kz); 
    Kstatic   = reshape(Kstatic,[param.Np, param.Nseg*param.Nshot, 1, 3]);
    
  % Compute uniform density compensation weights    
    param.filterFactor = 1;
    [wcut, wnocut] = dcf3Dradial(param.Np,param.Nseg*param.Nshot,param.filterFactor);

  % Apply same weights for all readout projections
    wcut   = repmat(  wcut(:), [1, param.Nseg * param.Nshot]);   
    wnocut = repmat(wnocut(:), [1, param.Nseg * param.Nshot]);   

    if param.dcf_cutFlag
        Wstatic = wcut;
    else
        Wstatic = wnocut;
    end
    
%--------------------------------------------------------------------------
% Here: your Fourier transform / reconstruction code using rawData, Kstatic
% and Wstatic
%--------------------------------------------------------------------------  

%% GROG 3D 

kdata = permute(rawData,[1 3 4 2]);
kx = reshape(kx,size(kx,1),size(kx,2)*size(kx,3));
ky = reshape(ky,size(kx,1),size(kx,2)*size(kx,3));
kz = reshape(kz,size(kx,1),size(kx,2)*size(kx,3));

s = 100;
kxs = kx(:,1:s);
kys = ky(:,1:s);
kzs = kz(:,1:s);
kdatas = kdata(:,1:s,:,:);

%%
% cc
clear Gx Gy Gz kdatas_tmp;

i = 1; step = 2;
for c = 1:step:30
    kdatas_tmp(:,:,:,i) = sum((kdatas(:,:,:,c:c+step-1)),4);
    i = i+1;
end

% figure; plot3(kxs,kys,kzs)
% [Gx,Gy,Gz] = GROG.get_Gx_Gy_Gz(kdatas,kxs,kys,kzs);
[Gx,Gy,Gz] = GROG.get_Gx_Gy_Gz(kdatas(:,:,:,[4:7,18:21,30]),kxs,kys,kzs);
% [Gx,Gy] = GROG.get_Gx_Gy(kdatas,kxs+i*kys);

% close all;
figure('Color','White','Position',[300 300 700 300]);
subplot(1,3,1)
imagesc(abs(Gx)); title('Gx'); 
subplot(1,3,2)
imagesc(abs(Gy)); title('Gy')
subplot(1,3,3)
imagesc(abs(Gz)); title('Gz')

%%
% G = GROG.init_3D(kdata,kx,ky,kz,Gx,Gy,Gz,0);
G = GROG.init_3D(kdatas(:,:,:,[4:7,18:21,30]),kxs,kys,kzs,Gx,Gy,Gz,0);

kref = GROG.interp_3D(kdatas(:,:,:,[4:7,18:21,30]),G,1);       % gridding kspace: radial -> cartesian

ref = squeeze(ifft3c_mri(kref));     % i2k linear recon

figure; 
for c=1:size(ref,4) 
    subplot(3,10,c); 
    imagesc(abs(squeeze(kref(:,:,112,c)))); 
end



%% CS recon
b1 = adapt_array_2d(ref); clear ref; % coil sensitivity map
b1 = single(b1/max(abs(b1(:))));

% plot coil sensitivity maps
figure('Color','White','Position',[300 300 700 300]);
for coil = 1:size(b1,3)
    subplot(2,4,coil); imagesc(abs(b1(:,:,coil))); axis off; colormap gray; 
    title(['Coil ', num2str(coil)])
end

% calculate the DCF for nyquist sampling - or 1/2 sub-Nyquist?
% R_slice = Nqu / #spokes_per_slice
Nqu = floor(bas/2*pi); % sampling rate at which full Nyquist is achieved? 
% use one set (time frame) of Nqu points
G = GROG.init(kdata(:,end-Nqu+1:end,:,:),Traj(:,end-Nqu+1:end),Gx,Gy,0);
% DCF is a ram-lak filter see doi:10.1002/mrm.27030.
DCF = reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1))]);
DCF = CropImg(DCF,nx,nx);

% sort the data
nt = floor(ntviews/nline);
Traj = reshape(Traj(:,1:nt*nline),[nx,nline,nt]);
kdata = reshape(kdata(:,1:nt*nline,:,:),[nx,nline,nt,nc]);
[nx,ntviews,nt,nc] = size(kdata);

% calculat weighting for iteration
G = GROG.init(kdata,Traj,Gx,Gy,0);
DCF_U = reshape(G.weight,[sqrt(size(G.weight,1)),sqrt(size(G.weight,1)),nt]);
DCF_U = CropImg(DCF_U,nx,nx);
DCF = repmat(DCF,[1,1,nt,nc]);
DCF_U = repmat(DCF_U,[1,1,1,nc]);
Weighting = DCF./DCF_U;

% grog
kdata = GROG.interp(kdata,G,1);
mask = single(kdata~=0);

%%
% plot data mask
figure('Color','White','Position',[300 300 600 500]);
for i = 1:20
    imagesc(abs(mask(:,:,i,1))); colormap gray; hold on; 
    pause(0.1)
end

param.y = kdata.*sqrt(Weighting);
param.E = Emat_GROG2D(mask,b1,Weighting);
clear mask Weighting b1

recon_cs=param.E'*param.y;
data=abs(single(gather(recon_cs/max(recon_cs(:)))));

param.TV=TV_Temp;
Weight1=0.06;
param.TVWeight=max(abs(recon_cs(:)))*Weight1;
param.nite = 7;param.display = 1;
clc
tic
for n=1:3
    recon_cs = CSL1NlCg(recon_cs,param);
end
time=toc/60
data(:,:,:,2)=abs(gather(single(recon_cs/max(recon_cs(:)))));

data=CropImg(data,bas,bas);

%% plot results

close all;
figure('Color','White','Position',[300 300 700 300]);
for n=1:size(data,3)
    subplot(1,2,1)
    imagesc(squeeze(data(:,:,n,1))); hold on; axis off;
    subplot(1,2,2)
    imagesc(squeeze(data(:,:,n,2))); hold on; axis off;
    pause(0.1);
end

    