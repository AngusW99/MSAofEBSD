%% Matlab script for factor/cluster analysis of EBSD data
%                                  - Angus J Wilkinson May 2018

% clear the decks by removing all variables and shutting all matlab figures
clear variables; close all;

%start a timer - so we can tell how long different tasks took
t1=clock;

%%  <<<< USER INPUTS >>>>
%   give location of dataset and specify pattern size to work on
%   array size can get very large so may need to resize to bin down from
%   experimental conditions


% small ferritic steel sample - from Merlin- from DMC
BaseFolder='F:\Ben\';
HDF5_filename='Demo_Ben.h5';
NumPCA=36; %number of factors to use
PCAGrid=[6 6]; %grid for plotting

PCA_radon=0; %do a PCA on the Radon
plot_pca=1; %plot the PCA

%%
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

%%
%load AstroEBSD
cur_folder=cd;
run('C:\Users\tbritton\Documents\GitHub\AstroEBSD\start_AstroEBSD');
cd(cur_folder);

%% Read in the data
pTime('Reading in h5 data',t1);
% open HDF5 file and read in some data on EBSD map
InputUser.HDF5_folder=BaseFolder;
InputUser.HDF5_file=HDF5_filename;
[ MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
%read the map & convert to area data
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

% this bit finds the extent of the cropped region of the SEM image that was actually mapped with EBSD
% this is not used properly - will need some cropped map data to test/fix
% the code for this!

MapSize1=max(MapData.YBeam)-min(MapData.YBeam)+1; % number of rows
MapSize2=max(MapData.XBeam)-min(MapData.XBeam)+1; % number of columns
SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));
figure;
subplot(1,2,1); imagesc(SEMImage); axis image off; colormap('gray'); title('SEM Image');
subplot(1,2,2); imagesc(transpose(reshape(MapData.RadonQuality,MapSize2,MapSize1)));axis image off;colormap('gray'); title('IQ');

%% Read in pattern data, correct, reshape
pTime('Single Pattern Checking',t1);
%load the first pattern to set up an array of patterns
[ RefPat ] = bReadEBSP(EBSPData,1);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
[R_EBSP0,R_RHO0]=radon(RefPat_cor*0+1,0:180);
[R_EBSP,R_RHO]=radon(RefPat_cor,0:180);
R_EBSP=R_EBSP./R_EBSP0;
ylims=(-size(RefPat_cor,1)/2:size(RefPat_cor,1)/2)+(size(R_EBSP,1)-1)/2;
R_EBSPr=R_EBSP(ylims,:);
PatSizeR=size(R_EBSPr,1);
PatSizeT=size(R_EBSPr,2);

%plot an example pattern
figure;
subplot(1,2,1); imagesc(RefPat); axis image; axis xy; axis tight; colormap('gray'); title('Input pattern');
subplot(1,2,2); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');

%% load the data
pTime('Loading Patterns',t1);
PatSizeW=size(RefPat_cor,1);
PatSizeH=size(RefPat_cor,2);
testArray=zeros(PatSizeW*PatSizeH,Data_InputMap.ypts,Data_InputMap.xpts);

if PCA_radon == 1
    testRadon=zeros(PatSizeR*PatSizeT,Data_InputMap.ypts,Data_InputMap.xpts);
end

for xi=1:Data_InputMap.xpts  %MapSize1*MapSize2
    parfor yi=1:Data_InputMap.ypts
        pattern_number=Data_InputMap.PMap(yi,xi);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        testArray(:,yi,xi)=reshape(RefPat_cor,PatSizeW*PatSizeH,1);
        
        if PCA_radon == 1
            [R_EBSP0,R_RHO0]=radon(RefPat_cor*0+1,0:180);
            [R_EBSP,R_RHO]=radon(RefPat_cor,0:180);
            R_EBSP=R_EBSP./R_EBSP0;
            ylims=(-size(RefPat_cor,1)/2:size(RefPat_cor,1)/2)+(size(R_EBSP,1)-1)/2;
            R_EBSPr=R_EBSP(ylims,:);
            testRadon(:,yi,xi)=reshape(R_EBSPr,PatSizeT*PatSizeR,1);
        end
    end
    pTime(['Completed row = ' int2str(xi)],t1);
end

%%
pTime('Reshaping Data',t1);
%reshape
testArray2=reshape(testArray,size(testArray,1),size(testArray,2)*size(testArray,3));

%flatten each pattern properly - this is important to get right
testmean=mean(testArray2,1); %calc the mean of each pattern
teststd=std(testArray2,0,1); %calc the std of each pattern
%flatten pattern by pattern
testArray2=(testArray2-repmat(testmean,[size(testArray2,1),1]))./repmat(teststd,[size(testArray2,1),1]);

%% PCA radon
%coeff = maps
%scores = patterns

if PCA_radon == 1
    pTime('PCA Radon',t1);
    testRadon2=reshape(testRadon,size(testRadon,1),size(testRadon,2)*size(testRadon,3));
    [coeffr,scorer,latent,~,explained]=pca(testRadon2, 'Centered',false, 'NumComponents',NumPCA);
    %rotate the factors
    [coeffVMr, RotVMr] = rotatefactors(coeffr(:,1:NumPCA)+mean(mean(coeffr)),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
    coeffVMr=coeffVMr-mean(mean(coeffr));
    scoreVMr=scorer(:,1:NumPCA)*RotVMr;
end

%% plot PCA on Radon
%plot the loading maps

if plot_pca == 1 && PCA_radon == 1
    figure;
    for n=1:NumPCA
        subplot(2,5,n);
        imagesc(reshape(coeffr(:,n),[PatSizeR,PatSizeT])); axis image; axis xy; axis tight; colormap('gray'); axis off;
    end
    figure;
    for n=1:NumPCA
        subplot(2,5,n);
        imagesc(reshape(scorer(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;
    end
    
    figure;
    for n=1:NumPCA
        subplot(2,5,n);
        imagesc(reshape(coeffVMr(:,n),[PatSizeR,PatSizeT])); axis image; axis xy; axis tight; colormap('gray'); axis off;
    end
    figure;
    for n=1:NumPCA
        subplot(2,5,n);
        imagesc(reshape(scoreVMr(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;
    end
end
%%  run PCA analysis to split data into basis patterns (eigenvectors) and signal strength in map (loadings)
%   this uses a user input choice on number of factors to retain (NumPCA)
pTime('PCA normal',t1);
[coeff,score,latent,~,explained]=pca(testArray2, 'Centered',false, 'NumComponents',NumPCA);
pTime('Rotating Factors',t1);
%rotate the factors
[coeffVM, RotVM] = rotatefactors(coeff(:,1:NumPCA),'Normalize','on','Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
scoreVM=score(:,1:NumPCA)*RotVM;

% sort the factor map
PCA_VM_Pat_n=zeros(PatSizeW,PatSizeH,NumPCA);
PCA_VM_Map_n=zeros([Data_InputMap.ypts Data_InputMap.xpts NumPCA]);
PCA_VM_num=zeros([Data_InputMap.ypts Data_InputMap.xpts]);
PCA_VM_val=zeros([Data_InputMap.ypts Data_InputMap.xpts]);

%sort the loadings
for n=1:NumPCA
    PCA_VM_Pat_n(:,:,n)=reshape(scoreVM(:,n),[PatSizeW,PatSizeH]);
    Map_n=reshape(coeffVM(:,n),[Data_InputMap.ypts Data_InputMap.xpts]);
    PCA_VM_Map_n(:,:,n)=Map_n;
    
    PCA_VM_num(abs(Map_n) > abs(PCA_VM_val)) = n;
    PCA_VM_val(abs(Map_n) > abs(PCA_VM_val)) = Map_n(abs(Map_n) > abs(PCA_VM_val));
end

%flip the loadings if needed
for n=1:NumPCA
    meanV=mean(PCA_VM_val(PCA_VM_num==n));
    if meanV < 0
        PCA_VM_Pat_n(:,:,n)=-PCA_VM_Pat_n(:,:,n);
        PCA_VM_Map_n(:,:,n)=-PCA_VM_Map_n(:,:,n);
        PCA_VM_val(PCA_VM_num==n)=-PCA_VM_val(PCA_VM_num==n);
    end
end
%% Plot the PCAs
if plot_pca == 1
    figure;
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(reshape(score(:,n),[PatSizeW,PatSizeH])); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
    end
    figure
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(reshape(coeff(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;title(int2str(n));
    end
    
    %plot these new factors
    figure;
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(PCA_VM_Pat_n(:,:,n)); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
    end
    figure
    for n=1:NumPCA
        subplot(PCAGrid(1),PCAGrid(2),n);
        imagesc(PCA_VM_Map_n(:,:,n)); axis image; axis ij; axis tight; axis off;title(int2str(n));
    end
    
    %plot the factor map
    figure;
    imagesc(PCA_VM_num); axis image; colormap('colorcube');
end
%% single point selection for pattern inspection
%{

%}

%% save the data
s=clock;
save_file=fullfile(BaseFolder,HDF5_filename(1:end-3),[int2str(s(1)) '_' int2str(s(2)) '_' int2str(s(3)) '_' int2str(s(4)) int2str(s(5)) int2str(s(6)) '.mat']);
if exist(fullfile(BaseFolder,HDF5_filename(1:end-3))) ~= 7
    mkdir(fullfile(BaseFolder,HDF5_filename(1:end-3)));
end
save(save_file);