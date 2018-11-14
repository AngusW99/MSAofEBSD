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
BaseFolder='G:\';
HDF5_filename='Si_100x_50um.h5';

NumPCA=8; %number of factors to use
PCAGrid=[4 2]; %grid for plotting

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
Settings_Cor.SquareCrop = 1; %make square the EBSP
%%
%load AstroEBSD
cur_folder=cd;
run('C:\Users\tbritton\Documents\GitHub\AstroEBSD\start_AstroEBSD');
run('C:\Users\tbritton\Documents\GitHub\rTemplateMatch\start_RTI');
cd(cur_folder);

%% Read in the data
pTime('Reading in h5 data',t1);
% open HDF5 file and read in some data on EBSD map
InputUser.Phase_Input  = {'Si'}; %Si, Ferrite

InputUser.HDF5_folder=BaseFolder;
InputUser.HDF5_file=HDF5_filename;
[ MapData,MicroscopeData,PhaseData,EBSD_DataInfo ]=bReadHDF5( InputUser );
%read the map & convert to area data
[Data_InputMap_Start] = EBSD_Map(MapData,MicroscopeData);
[ Data_InputMap ] = PC_square( EBSD_DataInfo, Data_InputMap_Start,Settings_Cor );
BinFile='C:\Users\tbritton\Documents\GitHub\rTemplateMatch\masterpatterns\Si_1024.bin'; isHex=0;
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
[ RefPat ] = bReadEBSP(EBSD_DataInfo,1);
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
        [ RefPat ] = bReadEBSP(EBSD_DataInfo,pattern_number);
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
%  save(save_file);

%% index the patterns


[screen_int,facedata] = Cube_Generate(BinFile,isHex);

%% Input variables
screensize=128; %in pixels - will force everything to square

%input the crystal symmentry you are using, cubic, hexagonal etc.
symmetry = 'cubic'; %'cubic'

%Set the number of iterations to do in the refinement step
iterations = 4;

%LPT size used
LPTsize = 500;
XCF_type=2; % 2 = FFT, NPD doesn't work properly...

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
Settings_CorA=Settings_Cor;
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
    EBSP_raw = bReadEBSP(EBSD_DataInfo,pnum);
    
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

%%
[ RefPat ] = bReadEBSP(EBSD_DataInfo,1);
[ RefPat_cor,PatternInfo] = EBSP_BGCor( RefPat,Settings_Cor );
PC_MSA=[PCfun_x(1),PCfun_y(1),PCfun_z(1)];

[ EBSP_av ] = EBSP_Gnom(PatternInfo,PC_MSA);
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( screensize, LPTsize );

Sampling_Freq=7;
cs = crystalSymmetry(symmetry,'mineral','ferrite');
[ library_G ] = SO3_rotmat_gen( cs,Sampling_Freq);
% generate the templatest
[ template_library ] = Library_Gen(EBSP_av,screen_int,isHex,library_G,XCF_type,SettingsXCF);
%send to the pool
template_library_par = parallel.pool.Constant(template_library);
screen_int_par=parallel.pool.Constant(screen_int);


%%
[ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_MSA);
parfor n=1:NumPCA

    [rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, PCA_VM_Pat_n(:,:,n), library_G, template_library_par.Value);
    % [rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, EBSP_cor, library_G, template_library);
    [rotmat_best,regout] = refine2(PCA_VM_Pat_n(:,:,n),PC_MSA,PC_MSA,EBSP_pat,rotmat_1st,Rx,Ry,Rz,shift_xFFT,shift_yFFT,SettingsXCF,SettingsXCF2,correction,screen_int_par.Value,isHex,LPTsize,iterations);
    %convert to euler angles
    eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
    MSA_Ang(:,n)=[eang_template regout]';
end

%% populate the EA data from the MSA fitting


for pat=1:Data_InputMap.max_pats
    Data_MSA_phi1(pat)=MSA_Ang(1,PCA_VM_num(pat));
    Data_MSA_phi2(pat)=MSA_Ang(3,PCA_VM_num(pat));
    Data_MSA_phi(pat)=MSA_Ang(2,PCA_VM_num(pat));
    props.x(pat)=Data_InputMap.XSample(pat);
    props.y(pat)=Data_InputMap.YSample(pat);
end

rot_msa=rotation('Euler',Data_MSA_phi1(:),Data_MSA_phi(:),Data_MSA_phi2(:),'ZXZ');
phases=ones(size(Data_MSA_phi1(:),1),1);
EBSD_mas=EBSD(rot_msa,phases,cs,'options',props);

setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');
phase='ferrite';
oM=ipfTSLKey(EBSD_mas(phase));

figure;
subplot(1,3,1)
oM.inversePoleFigureDirection = xvector;
plot(EBSD_mas(phase), oM.orientation2color(EBSD_mas(phase).orientations),'parent',gca)
title('Astro - X')

subplot(1,3,2)
oM.inversePoleFigureDirection = yvector;
plot(EBSD_mas(phase), oM.orientation2color(EBSD_mas(phase).orientations),'parent',gca)
title('Astro - Y')

subplot(1,3,3)
oM.inversePoleFigureDirection = zvector;
plot(EBSD_mas(phase), oM.orientation2color(EBSD_mas(phase).orientations),'parent',gca)
title('Astro - Z')

%% Refine the OriSolution
AstroRefined=zeros(7,Data_InputMap.max_pats);
tic
parfor pat=1:Data_InputMap.max_pats
%     pnum=Data_InputMap.PMap(pat);
%     EBSP_raw = bReadEBSP(EBSD_DataInfo,pnum);
%     %build the geometry
%     [ PatternCor,PatternInfo ] = EBSP_BGCor( EBSP_raw,Settings_CorX );
    PatternCor=reshape(testArray2(:,pat),128,128);
    
    PC_X=PCfun_x(1)+Data_InputMap.XBeam_Map(pat)*PCfun_x(2)+Data_InputMap.XBeam_Map(pat)*PCfun_x(3);
    PC_Y=PCfun_y(1)+Data_InputMap.XBeam_Map(pat)*PCfun_y(2)+Data_InputMap.YBeam_Map(pat)*PCfun_y(3);
    PC_Z=PCfun_z(1)+Data_InputMap.XBeam_Map(pat)*PCfun_z(2)+Data_InputMap.YBeam_Map(pat)*PCfun_z(3);
    PC_pat=[PC_X,PC_Y,PC_Z];
    
    [ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);
    RotMat_Astro=conv_EA_to_G([Data_MSA_phi1(pat),Data_MSA_phi(pat),Data_MSA_phi2(pat)])*Detector_tilt;
    [rotmat_astro_refined,regout] = refine2(PatternCor,PC_pat,PC_pat,EBSP_pat,RotMat_Astro,Rx,Ry,Rz,0,0,SettingsXCF,SettingsXCF2,correction,screen_int_par.Value,isHex,LPTsize,iterations+1);
    eang_astro_refined=conv_G_to_EA(rotmat_astro_refined*inv(Detector_tilt));
    
    AstroRefined(:,pat)=[eang_astro_refined regout]';
end
toc
%%
rot_astroref=rotation('Euler',AstroRefined(1,:),AstroRefined(2,:),AstroRefined(3,:),'ZXZ');
EBSD_astroref=EBSD(rot_astroref,phases,cs,'options',props);

figure;
subplot(1,3,1)
oM.inversePoleFigureDirection = xvector;
plot(EBSD_astroref(phase), oM.orientation2color(EBSD_astroref(phase).orientations),'parent',gca)
title('Astro - X')

subplot(1,3,2)
oM.inversePoleFigureDirection = yvector;
plot(EBSD_astroref(phase), oM.orientation2color(EBSD_astroref(phase).orientations),'parent',gca)
title('Astro - Y')

subplot(1,3,3)
oM.inversePoleFigureDirection = zvector;
plot(EBSD_astroref(phase), oM.orientation2color(EBSD_astroref(phase).orientations),'parent',gca)
title('Astro - Z')

%% clean up
gbThreshold=5*degree;
[grains,ebsd.grainId] = calcGrains(EBSD_astroref('indexed'),'angle',gbThreshold);

smallgrains=grains(grains.grainSize<=4);
sg_x=smallgrains.x;
sg_y=smallgrains.y;
n_sg=size(sg_y,1);

AstroRefined2=AstroRefined;

for n=1:n_sg
    %find position
    r_data=sqrt((props.x-sg_x(n)).^2+(props.y-sg_y(n)).^2);
    ix=find(r_data<2*Data_InputMap.XStep);
    [~,pat]=min(r_data(:));
    
    PatternCor=reshape(testArray2(:,pat),128,128);
    
    PC_X=PCfun_x(1)+Data_InputMap.XBeam_Map(pat)*PCfun_x(2)+Data_InputMap.XBeam_Map(pat)*PCfun_x(3);
    PC_Y=PCfun_y(1)+Data_InputMap.XBeam_Map(pat)*PCfun_y(2)+Data_InputMap.YBeam_Map(pat)*PCfun_y(3);
    PC_Z=PCfun_z(1)+Data_InputMap.XBeam_Map(pat)*PCfun_z(2)+Data_InputMap.YBeam_Map(pat)*PCfun_z(3);
    PC_pat=[PC_X,PC_Y,PC_Z];
    
    [ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);

    
    regout=zeros(size(ix,2),4);
    eang_astro_refix=zeros(size(ix,2),3);
    parfor p=1:size(ix,2)
        p_try=ix(p);
        RotMat_Astro=conv_EA_to_G([Data_MSA_phi1(p_try),Data_MSA_phi(p_try),Data_MSA_phi2(p_try)])*Detector_tilt;
        [rotmat_astro_refined,regout(p,:)] = refine2(PatternCor,PC_pat,PC_pat,EBSP_pat,RotMat_Astro,Rx,Ry,Rz,0,0,SettingsXCF,SettingsXCF2,correction,screen_int_par.Value,isHex,LPTsize,iterations+1);
        eang_astro_refix(p,:)=conv_G_to_EA(rotmat_astro_refined*inv(Detector_tilt));
    end
    [~,ixp]=max(regout(:,4));
    AstroRefined2(:,pat)=[eang_astro_refix(ixp,:) regout(ixp,:)]';
end
%%
rot_astroref2=rotation('Euler',AstroRefined2(1,:),AstroRefined2(2,:),AstroRefined2(3,:),'ZXZ');
EBSD_astroref2=EBSD(rot_astroref2,phases,cs,'options',props);

figure;
subplot(1,3,1)
oM.inversePoleFigureDirection = xvector;
plot(EBSD_astroref2(phase), oM.orientation2color(EBSD_astroref2(phase).orientations),'parent',gca)
title('RAstro - X')

subplot(1,3,2)
oM.inversePoleFigureDirection = yvector;
plot(EBSD_astroref2(phase), oM.orientation2color(EBSD_astroref2(phase).orientations),'parent',gca)
title('RAstro - Y')

subplot(1,3,3)
oM.inversePoleFigureDirection = zvector;
plot(EBSD_astroref2(phase), oM.orientation2color(EBSD_astroref2(phase).orientations),'parent',gca)
title('RAstro - Z')