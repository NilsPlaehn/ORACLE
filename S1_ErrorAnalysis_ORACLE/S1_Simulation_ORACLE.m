
% Description: Montecarlo simulation for error characteristic of quantitative ORACLE in  presence of finite noise

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
% nvari:    number of variations of complex normal distributed noise (e.g. Monte Carlo runs)
% SNR:      mean signal per standard deviation of noise
clear all
TR      = 5/1000;  
TE      = TR/2;
nPC     = 20; 
theta   = pi/(2*nPC); 
alpha   = 15.*pi/180;
M0      = 1;

nT1     = 105;
nT2     = 100; 

T1vec      = linspace(20,3000,nT1)./1000;
T2vec      = linspace(1,1000,nT2)./1000;

phit    = linspace(0,2*pi,nPC+1);
phi     = phit(1:nPC);

nvari = 10^5;
SNR   = 100;  

%% 1.2) Determine mean signal for ROI
bSSFPSignal = zeros(nPC,nT2,nT1);
for indT2 = 1:nT2
    for indT1 = 1:nT1
        T1 = T1vec(indT1);
        T2 = T2vec(indT2);
        for indPC = 1:nPC
            bSSFPSignal(indPC,indT2,indT1) = S1_bSSFP_Profile_Generation(M0,T1,T2,alpha,phi(indPC),TR,TE,theta);%bSSFPSignal_generator_fct(T1w,T2w,offw,TR,alpha,nPC);
        end
    end
end

Amax         = mean(mean(mean(abs(bSSFPSignal(:,:,:)))));
%% Determine T1 and T2 error estimation in presence of noise 

T1error = zeros(nT1,nT2);
T2error = zeros(nT1,nT2);

M0error = zeros(nT1,nT2);
B0error = zeros(nT1,nT2);

for indT1 = 1:nT1
    for indT2 = 1:nT2
        T1 = T1vec(indT1);
        T2 = T2vec(indT2);
        
        T1tmp  = zeros(1,nvari);
        T2tmp  = zeros(1,nvari);
        T1tmp2 = zeros(1,nvari);
        T2tmp2 = zeros(1,nvari);

        M0tmp = zeros(1,nvari);
        M0tmp2 = zeros(1,nvari);

        B0tmp  = zeros(1,nvari);
        B0tmp2 = zeros(1,nvari);
        for inditer = 1:nvari 
            profile0 = zeros(1,nPC);
            theta =1;  %rand(1)*100;
            for indPC = 1:nPC
                profile0(indPC) = S1_bSSFP_Profile_Generation(M0,T1,T2,alpha,phi(indPC),TR,TE,theta);
            end
            % sqrt(2) because abs(1+1i)=sqrt(2) 
            profile           = profile0+Amax/SNR/sqrt(2)*(randn(size(profile0))+1i*randn(size(profile0)));

            [T1_est,T2_est, B0_est, M0_est] =S1_ORACLE_fct(profile,TR,alpha);%
            T1tmp(inditer) = T1_est/T1;
            T2tmp(inditer) = T2_est/T2;
            M0tmp(inditer) = M0_est/M0;
            B0tmp(inditer) = angle(exp(1i*(B0_est)));
        end

        
        T1error(indT1,indT2) = std(T1tmp)*100;
        T2error(indT1,indT2) = std(T2tmp)*100;
        M0error(indT1,indT2) = std(M0tmp)*100;
        B0error(indT1,indT2) = std(B0tmp);

    end
    indT1
end



%% sff

close all


%% Saf
Plot_T1(2,T2vec,T1vec,T1error)
Plot_T2(4,T2vec,T1vec,T2error)
Plot_M0(5,T2vec,T1vec,M0error)

%% SAF
Plot_B0(7,T2vec,T1vec,B0error(:,:))


%% FUNCTIONS ________________________________________________
%%           ________________________________________________
%%           ________________________________________________
%%           ________________________________________________





function Plot_M0(indFigure,T2,T1,M0error)
    levels =0:0.1:3;
    f = figure(indFigure);
    c = contourf(T2,T1,M0error,levels,'edgecolor','none');
    h = colorbar;%('northoutside');
   % colormap(inferno)
   % set(gca,'ColorScale','log');
    ylabel('T_1 in s','FontSize',20)
    xlabel('T_2 in s','FontSize',20)
    title('|\Delta M_0| in %')
    ylabel(h,'relative error of M_0 in %','FontSize',20)
    h.Ticks = levels; 
    h.TickLength = 0.068;
  
    
    ax = gca; 
    
    ax.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];
    
    ax.YAxis.MinorTickValues =  (100:100:3000)/1000;
    ax.YAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.XAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
    ax.XAxis.MinorTickValues = (100:100:1000)/1000;
    ax.XAxis.TickValues = [200 400 600 800 1000]/1000;
    % 
    ax.ZMinorTick = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.7; % maximum line opacity
    set(ax,'FontSize',16)
    ax.MinorGridColor = 'k';%[.7 .7 .7];%'w';%[0.4660 0.6740 0.1880];
    ax.MinorGridAlpha = 0.25;
    
    f.Position = [1000 300 550 420];
end

function Plot_B0(indFigure,T2,T1,B0error)
    levels = 0.01:0.1:5;
    f = figure(indFigure);
    c = contourf(T2,T1,180/pi*B0error,levels,'edgecolor','none');
    h = colorbar;%('northoutside');
  %  colormap(inferno)
    clim([0.01 5])
    set(gca,'ColorScale','log');
    ylabel('T_1 in s','FontSize',20)
    xlabel('T_2 in ms','FontSize',20)
    %title('|\Delta T_2| in %')
    title('|\Delta \phi_1| in ()°')
    ylabel(h,'error of \phi_1 in ()°','FontSize',20)
    h.Ticks = [0.1 0.2 0.4 1 2 3];
    h.TickLength = 0.068;
    
    
    %h.Color = [0.95 0.95 0.95];
    
    ax = gca; 
    
    ax.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];
    
    ax.YAxis.MinorTickValues =  (100:100:3000)/1000;
    ax.YAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.XAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
    ax.XAxis.MinorTickValues = (100:100:1000)/1000;
    ax.XAxis.TickValues = [200 400 600 800 1000]/1000;
    % 
    ax.ZMinorTick = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.7; % maximum line opacity
    set(ax,'FontSize',16)
    ax.MinorGridColor = 'k';%[.7 .7 .7];%'w';%[0.4660 0.6740 0.1880];
    ax.MinorGridAlpha = 0.25;
    
    f.Position = [200 300 550 420];
    ax
end

function Plot_T1(indFigure,T2,T1,Tierror)

    levels = 0:2:14;
    f = figure(indFigure);
    c = contourf(T2,T1,Tierror,levels,'edgecolor','none');
    h = colorbar;%('northoutside');
    %colormap(hot)
    %set(gca,'ColorScale','log');
    ylabel('T_1 in s','FontSize',20)
    xlabel('T_2 in s','FontSize',20)
    %set(gca,'YScale','log')
    % set(gca,'XScale','log')
    title('|\Delta T_1| in %')
    ylabel(h,'relative error of T_1 in %','FontSize',20)
    h.Ticks = levels; 
    h.TickLength = 0.068;
    clim([0 14])
    
    %h.Color = [0.95 0.95 0.95];
    
    ax = gca; 
    
    ax.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];
    
    ax.YAxis.MinorTickValues =  (100:100:3000)/1000;
    ax.YAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.XAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
    ax.XAxis.MinorTickValues = (100:100:1000)/1000;
    ax.XAxis.TickValues = [200 400 600 800 1000]/1000;
    % 
    ax.ZMinorTick = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.7; % maximum line opacity
    set(ax,'FontSize',16)
    ax.MinorGridColor = 'k';%[.7 .7 .7];%'w';%[0.4660 0.6740 0.1880];
    ax.MinorGridAlpha = 0.25;
    
    f.Position = [1000 300 550 420];
end

function Plot_T2(indFigure,T2,T1,Tierror)

    levels = 0:2:14;
    f = figure(indFigure);
    c = contourf(T2,T1,Tierror,levels,'edgecolor','none');
    h = colorbar;%('northoutside');
    %colormap(hot)
    %set(gca,'ColorScale','log');
    ylabel('T_1 in s','FontSize',20)
    xlabel('T_2 in s','FontSize',20)
    %set(gca,'YScale','log')
    % set(gca,'XScale','log')
    title('|\Delta T_2| in %')
    ylabel(h,'relative error of T_2 in %','FontSize',20)
    h.Ticks = levels; 
    h.TickLength = 0.068;
    clim([0 14])
    
    %h.Color = [0.95 0.95 0.95];
    
    ax = gca; 
    
    ax.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];
    
    ax.YAxis.MinorTickValues =  (100:100:3000)/1000;
    ax.YAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.XAxis.TickLength = [0.0100 0.0250]*1.5;
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
    ax.XAxis.MinorTickValues = (100:100:1000)/1000;
    ax.XAxis.TickValues = [200 400 600 800 1000]/1000;
    % 
    ax.ZMinorTick = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.GridLineStyle = ':';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.7; % maximum line opacity
    set(ax,'FontSize',16)
    ax.MinorGridColor = 'k';%[.7 .7 .7];%'w';%[0.4660 0.6740 0.1880];
    ax.MinorGridAlpha = 0.25;
    
    f.Position = [1000 300 550 420];
end