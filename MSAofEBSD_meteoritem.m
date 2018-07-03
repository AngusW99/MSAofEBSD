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
HDF5_filename='IronMeteorite_Gibeon_ArgusOptimisationPaper_Map1.h5';
NumPCA=10;  % user input of number of basis patterns

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
Settings_Cor.size=256; %image width
Settings_Cor.SplitBG=1; %split chip fix
%%
%load AstroEBSD
cur_folder=cd;
run('C:\Users\tbritton\Documents\GitHub\AstroEBSD\start_AstroEBSD');
cd(cur_folder);
pTime('Starting Read of Data',t1);
%%  Setting up to read in data
SaveFolder=[BaseFolder 'VMfromHDF5'];
mkdir(SaveFolder);
SaveFolder=[SaveFolder '\'];
mkdir([SaveFolder 'Patterns']);
mkdir([SaveFolder 'Maps']);
%%
pTime('Reading in h5 data',t1);
% open HDF5 file and read in some data on EBSD map
InputUser.HDF5_folder=BaseFolder;
InputUser.HDF5_file=HDF5_filename;
[ MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
%read the map & convert to area data
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);

% this bit finds the extent of the cropped region of the SEM image that was
% actually mapped with EBSD
MapSize1=max(MapData.YBeam)-min(MapData.YBeam)+1; % number of rows
MapSize2=max(MapData.XBeam)-min(MapData.XBeam)+1; % number of columns
Image=transpose(MicroscopeData.SEMImage(:,:,1));
figure;
subplot(1,2,1); imagesc(Image);axis image off;colormap('gray'); title('SEM Image');
subplot(1,2,2); imagesc(transpose(reshape(MapData.RadonQuality,MapSize2,MapSize1)));axis image off;colormap('gray'); title('IQ');


%% Read the EDX spectra
channum=2048;
pattern_number=1;
[EDSData_cor,EDSData_raw ] = bReadEDX(EBSPData,pattern_number,channum);

%% Read in pattern data, correct, reshape

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
subplot(1,2,1); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');

PatSizeW=size(RefPat_cor,1);
PatSizeH=size(RefPat_cor,2);

testArray=zeros(PatSizeW*PatSizeH,Data_InputMap.ypts,Data_InputMap.xpts);
testRadon=zeros(PatSizeR*PatSizeT,Data_InputMap.ypts,Data_InputMap.xpts);

for xi=1:Data_InputMap.xpts  %MapSize1*MapSize2
    parfor yi=1:Data_InputMap.ypts
        pattern_number=Data_InputMap.PMap(yi,xi);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        testArray(:,yi,xi)=reshape(RefPat_cor,PatSizeW*PatSizeH,1);
        
        [R_EBSP0,R_RHO0]=radon(RefPat_cor*0+1,0:180);
        [R_EBSP,R_RHO]=radon(RefPat_cor,0:180);
        R_EBSP=R_EBSP./R_EBSP0;
        ylims=(-size(RefPat_cor,1)/2:size(RefPat_cor,1)/2)+(size(R_EBSP,1)-1)/2;
        R_EBSPr=R_EBSP(ylims,:);
        
        testRadon(:,yi,xi)=reshape(R_EBSPr,PatSizeT*PatSizeR,1);
    end
    pTime(['Completed row = ' int2str(xi)],t1);
end

%%
%reshape
testArray2=transpose(reshape(testArray,size(testArray,1),size(testArray,2)*size(testArray,3)));

%flatfield by mean pattern
testArray2=testArray2-repmat(trimmean(testArray2,90),[size(testArray2,1),1]);

testRadon2=transpose(reshape(testRadon,size(testRadon,1),size(testRadon,2)*size(testRadon,3)));
testRadon2=testRadon2-repmat(trimmean(testRadon2,90),[size(testRadon2,1),1]);
%% PCA on Radon
[coeffr,scorer,latent,~,explained]=pca(testRadon2, 'Centered',false, 'NumComponents',NumPCA);

%plot the loading maps
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

%rotate the factors
[coeffVMr, RotVMr] = rotatefactors(coeffr(:,1:NumPCA)+mean(mean(coeffr)),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
coeffVMr=coeffVMr-mean(mean(coeffr));
scoreVMr=scorer(:,1:NumPCA)*RotVMr;

%plot maps
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

%%  run PCA analysis to split data into basis patterns (eigenvectors) and signal strength in map (loadings)
%   this uses a user input choice on number of factors to retain (NumPCA)
    
[coeff,score,latent,~,explained]=pca(testArray2, 'Centered',false, 'NumComponents',NumPCA);
pTime('PCA complete',t1);
%%
figure;
for n=1:NumPCA
    subplot(2,5,n);
    imagesc(reshape(coeff(:,n),[PatSizeW,PatSizeH])); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
end
figure;
for n=1:NumPCA
    subplot(2,5,n);
    imagesc(reshape(score(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;title(int2str(n));
end

%rotate the factors
[coeffVM, RotVM] = rotatefactors(coeff(:,1:NumPCA),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));

scoreVM=score(:,1:NumPCA)*RotVM;

%sort positive terms
for n=1:NumPCA
    if sum(sum(coeffVM(:,n))) < 0
        coeffVM(:,n)=-coeffVM(:,n);
        scoreVM(:,n)=-scoreVM(:,n);
    end
end
pTime('Factors rotated',t1);

%plot these new factors
figure;
for n=1:NumPCA
    subplot(2,5,n);
    imagesc(reshape(coeffVM(:,n),[PatSizeW,PatSizeH])); axis image; axis xy; axis tight; colormap('gray'); axis off; title(int2str(n));
end

%%
figure;
for n=1:NumPCA
    subplot(2,5,n);
    imagesc(reshape(scoreVM(:,n),[Data_InputMap.ypts Data_InputMap.xpts])); axis image; axis ij; axis tight; axis off;title(int2str(n)); caxis([-100 500]);
end
%%
%sort a map of the factors
scoremap_n=0*reshape(abs(scoreVM(:,1)),[Data_InputMap.ypts Data_InputMap.xpts]);
scoremap_v=0*reshape(abs(scoreVM(:,1)),[Data_InputMap.ypts Data_InputMap.xpts]);

% for n=1:NumPCA
%     scoremap=reshape(abs(scoreVM(:,n)),[Data_InputMap.ypts Data_InputMap.xpts]);
%     scoremap_n(scoremap>scoremap_v)=n;
%     scoremap_v(scoremap>scoremap_v)=scoremap(scoremap>scoremap_v);
% end

for pi=1:(Data_InputMap.ypts*Data_InputMap.xpts)
    [yi,xi]=ind2sub([Data_InputMap.ypts,Data_InputMap.xpts],pi);
    si=scoreVM(pi,:);
    si=abs(si);
    [mv,iv]=max(si);
    si(iv)=[];
    scoremap_v(yi,xi)=mv./std(abs(si));
    scoremap_n(yi,xi)=iv;
end


figure;

subplot(1,2,1); imagesc(scoremap_v); axis image; axis xy; axis ij; caxis([0 10])
subplot(1,2,2); imagesc(scoremap_n); axis image; axis xy; axis ij;

%%
%save the data
s=clock;
save_file=fullfile(BaseFolder,HDF5_filename(1:end-4),[int2str(s(1)) '_' int2str(s(2)) '_' int2str(s(3)) '_' int2str(s(4)) int2str(s(5)) int2str(s(6)) '.mat']);
save(save_file);