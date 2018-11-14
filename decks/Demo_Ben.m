%% Matlab script for factor/cluster analysis of EBSD data
%                                  Derived from work by - Angus J Wilkinson May 2018
%                                  Adapted - Ben Britton July 2018

% clear the decks by removing all variables and shutting all matlab figures
clear variables; close all;

%start a timer - so we can tell how long different tasks took
t1=clock;

%% File inputs
BaseFolder='F:\Ben\';
HDF5_filename='Demo_Ben.h5';
Astro_Loc='C:\Users\tbritton\Documents\GitHub\AstroEBSD\start_AstroEBSD';
NumPCA=36; %number of factors to use

PCA_radon=0; %do a PCA on the Radon

plot_pca=1; %plot the PCA
PCAGrid=[6 6]; %grid for plotting (needs to add up, or be greater than NumPCA)

%% AstroEBSD based background correction
%background correction
Settings_Cor.gfilt=1; %use a low pass filter (do you mean high pass?)
Settings_Cor.gfilt_s=4; %low pass filter sigma

%radius mask
Settings_Cor.radius=0; %use a radius mask
Settings_Cor.radius_frac=0.95; %fraction of the pattern width to use as the mask

%hot pixel
Settings_Cor.hotpixel=1; %hot pixel correction
Settings_Cor.hot_thresh=1000; %hot pixel threshold

%resize
Settings_Cor.resize=1; %resize correction
Settings_Cor.size=128; %image height
Settings_Cor.SplitBG=1; %split chip fix
Settings_Cor.Square=1; %square crop

%% Run the scripts

%load the MSA files from the deck folder
[MSA_Loc,~,~]=fileparts(cd);
addpath(MSA_Loc);
addpath(fullfile(MSA_Loc,'bin'));
addpath(cd);

%load AstroEBSD - needed for functions
cur_folder=cd;
run(Astro_Loc);
cd(cur_folder);

LoadMap; %load the mapdata
LoadEBSP; %load the EBSP data - this may take a while..
rPCA; %run the PCA analysis
pPCA; %plot the PCA analysis

%% save the data
s=clock;
save_file=fullfile(BaseFolder,HDF5_filename(1:end-3),[int2str(s(1)) '_' int2str(s(2)) '_' int2str(s(3)) '_' int2str(s(4)) int2str(s(5)) int2str(s(6)) '.mat']);
if exist(fullfile(BaseFolder,HDF5_filename(1:end-3))) ~= 7
    mkdir(fullfile(BaseFolder,HDF5_filename(1:end-3)));
end
save(save_file);

%% index
BinFile='C:\Users\tbritton\Documents\GitHub\rTemplateMatch\masterpatterns\Ferrite_1024.bin'; isHex=0;
InputUser.Phase_Input  = {'Ferrite'}; %Si, Ferrite

screensize=128; %in pixels - will force everything to square
iterations = 5; %as calibrated previously
LPTsize = 500; %LPT size to use
XCF_type=2; % 2 = FFT, NPD doesn't work properly...


%% load the H5 file
%Unpack the data
%update the pattern centre for square cropping
[ Data_InputMap2 ] = PC_square( EBSPData, Data_InputMap,Settings_Cor );

%% build the rotation matricies
%generate the rotation matrix 
%they are radians
Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];

%Define the tilt of the detector
Detector_tilt = Rx(MicroscopeData.TotalTilt);

%% Astro Index - for pattern centre
tic
Astro_FP= 'C:\Users\tbritton\Documents\GitHub\AstroEBSD'; %please change this to your astro file path
InputUser.Phase_Folder = fullfile(Astro_FP,'phases');

InputUser.PatternPhase = 1; %phase number for this pattern
InputUser.PCSearch = 1;     %find the pattern centre - for single patterns


%Peak Finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

%peak hunt
Settings_Rad.max_peaks=10; %max number of peaks to return
Settings_Rad.num_peak=15; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=4; %search size in theta steps
Settings_Rad.rho_search_per=0.05; %radon search in fractions - describes the allowable band width
Settings_Rad.min_peak_width=0.002; %seperation of the peak width, in pixels

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( InputUser.Phase_Input,InputUser.Phase_Folder );

%background correction
Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=5; %low pass filter sigma

%radius mask
Settings_CorX.radius=0; %use a radius mask
Settings_CorX.radius_frac=0.85; %fraction of the pattern width to use as the mask

%hold pixel
Settings_CorX.hotpixel=1; %hot pixel correction
Settings_CorX.hot_thresh=1000; %hot pixel threshold

%resize
Settings_CorX.resize=1; %resize correction
Settings_CorX.size=screensize; %image height

Settings_CorX.RealBG=0; %use a real BG
Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG

Settings_CorX.SquareCrop = 1; %make square the EBSP

Settings_CorX.SplitBG=1; %deal with a split screen

%fit the PC grid
PC_grid=[6 8]; %[y x];

%build the grid
[PC_xgrid,PC_ygrid]=meshgrid(1:PC_grid(2),1:PC_grid(1));
PC_xgrid=floor((Data_InputMap.xpts/(PC_grid(2)+1))*PC_xgrid);
PC_ygrid=floor((Data_InputMap.ypts/(PC_grid(1)+1))*PC_ygrid);

PC_xgrid=PC_xgrid(:);
PC_ygrid=PC_ygrid(:);

num_grid=size(PC_xgrid,1);

%settings
Settings_CorA=Settings_CorX;
Settings_CorA.radius=1;
Settings_CorA.size=512;

%set up the GA
num_P=1;

PC_GA_options = optimoptions('ga');
PC_GA_options.FunctionTolerance=1E-3;
PC_GA_options.UseParallel=0;
PC_GA_options.MaxGenerations=15;
PC_GA_options.PopulationSize=30;
PC_GA_options.MaxStallGenerations=20;
PC_GA_options.Display='off';
PC_found=zeros(3,num_grid);

parfor g=1:num_grid
    pnum=Data_InputMap.PMap(PC_ygrid(g),PC_xgrid(g));
    EBSP_raw = bReadEBSP(EBSPData,pnum);
    
    Settings_PCin_start=[Data_InputMap.PCX(pnum) Data_InputMap.PCY(pnum) Data_InputMap.DD(pnum)];
    Settings_PCin_range=[0.1 0.1 0.1];
    
    [ PatternCor,PatternInfo ] = EBSP_BGCor( EBSP_raw,Settings_CorA );
    % radon convert & Peak ID
    [ EBSP_Onep_Peak_Centre,EBSP_Onep_Single_Peak_Set_All,EBSP_Onep_Peak_Set_All,...
        EBSP_Onep_R_EBSP,EBSP_Onep_R_Edge,EBSP_Onep_R_rho,EBSP_Onep_R_theta ] ...
        = EBSP_RadHunt( PatternCor,Settings_Rad);
    % find the pattern centre
    
    PC_GA_ub=Settings_PCin_start+Settings_PCin_range;
    PC_GA_lb=Settings_PCin_start-Settings_PCin_range;
    
    EBSP_Onep_PC_out=zeros(3,Phase_Num);
    EBSP_Onep_PC_err=zeros(Phase_Num,1);
    
    FitFunc = @(PC_test) PC_GAOpt( PC_test,EBSP_Onep_Peak_Centre,PatternInfo.size,Crystal_LUT,Crystal_UCell,num_P);
    [EBSP_Onep_PC_out(:,num_P), EBSP_Onep_PC_err(num_P)] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
    
    PC_found(:,g)=EBSP_Onep_PC_out(:,num_P);
    
    
end

xyg=[PC_xgrid PC_ygrid];
PCfun_x=robustfit(xyg,PC_found(1,:));
PCfun_y=robustfit(xyg,PC_found(2,:));
PCfun_z=robustfit(xyg,PC_found(3,:));

%to eval
PC_X=PCfun_x(1)+PC_xgrid*PCfun_x(2)+PC_ygrid*PCfun_x(3);
PC_Y=PCfun_y(1)+PC_xgrid*PCfun_y(2)+PC_ygrid*PCfun_y(3);
PC_Z=PCfun_z(1)+PC_xgrid*PCfun_z(2)+PC_ygrid*PCfun_z(3);
disp('PC Model built');
toc

%% pattern match each master
PC_X_av=mean(PC_X(:));
PC_Y_av=mean(PC_Y(:));
PC_Z_av=mean(PC_Z(:));

