%% Matlab script for factor/cluster analysis of EBSD data
%                                  Derived from work by - Angus J Wilkinson May 2018
%                                  Adapted - Ben Britton July 2018

% clear the decks by removing all variables and shutting all matlab figures
clear variables; close all;

%start a timer - so we can tell how long different tasks took
t1=clock;

%% File inputs
BaseFolder='F:\Ben\';
HDF5_filename='IronMeteorite_Gibeon_ArgusOptimisationPaper_Map1.h5';
Astro_Loc='C:\Users\tbritton\Documents\GitHub\AstroEBSD\start_AstroEBSD';
NumPCA=200; %number of factors to use

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
Settings_Cor.size=64; %image height
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