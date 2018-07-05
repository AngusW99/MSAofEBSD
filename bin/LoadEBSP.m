%% Read in pattern data, correct, reshape
pTime('Single Pattern Checking',t1);

%load the first pattern to set up an array of patterns
[ RefPat ] = bReadEBSP(EBSPData,1);
%apply the correction
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );

if PCA_radon == 1
    %radon transform a blank image
    [R_EBSP0,R_RHO0]=radon(RefPat_cor*0+1,0:180);
    %radon transform the real image
    [R_EBSP,R_RHO]=radon(RefPat_cor,0:180);
    %normalise due to path length issues
    R_EBSP=R_EBSP./R_EBSP0;
    
    %crop the Radon
    ylims=(-size(RefPat_cor,1)/2:size(RefPat_cor,1)/2)+(size(R_EBSP,1)-1)/2;
    R_EBSPr=R_EBSP(ylims,:);
    PatSizeR=size(R_EBSPr,1);
    PatSizeT=size(R_EBSPr,2);
    %create the big radon array
    testRadon=zeros(PatSizeR*PatSizeT,Data_InputMap.ypts,Data_InputMap.xpts);
end

%% plot data
if PCA_radon == 0
    %plot an example pattern
    figure;
    subplot(1,2,1); imagesc(RefPat); axis image; axis xy; axis tight; colormap('gray'); title('Input pattern');
    subplot(1,2,2); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');
else
    %plot an example pattern
    figure;
    subplot(1,3,1); imagesc(RefPat); axis image; axis xy; axis tight; colormap('gray'); title('Input pattern');
    subplot(1,3,2); imagesc(RefPat_cor); axis image; axis xy; axis tight; colormap('gray'); title('Corrected pattern');
    subplot(1,3,3); imagesc(testRadon); axis image; axis xy; axis tight; colormap('gray'); title('Radon pattern');
end

%% load the data
pTime('Loading Patterns',t1);
PatSizeW=size(RefPat_cor,1);
PatSizeH=size(RefPat_cor,2);
testArray=zeros(PatSizeW*PatSizeH,Data_InputMap.ypts,Data_InputMap.xpts);

for xi=1:Data_InputMap.xpts  %MapSize1*MapSize2
    parfor yi=1:Data_InputMap.ypts
        pattern_number=Data_InputMap.PMap(yi,xi);
        [ RefPat ] = bReadEBSP(EBSPData,pattern_number);
        [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
        testArray(:,yi,xi)=reshape(RefPat_cor,PatSizeW*PatSizeH,1);
        
        if PCA_radon == 1
            [R_EBSP,R_RHO]=radon(RefPat_cor,0:180);
            R_EBSP=R_EBSP./R_EBSP0;
            ylims=(-size(RefPat_cor,1)/2:size(RefPat_cor,1)/2)+(size(R_EBSP,1)-1)/2;
            R_EBSPr=R_EBSP(ylims,:);
            testRadon(:,yi,xi)=reshape(R_EBSPr,PatSizeT*PatSizeR,1);
        end
    end
    pTime(['Completed row = ' int2str(xi)],t1);
end

%% Reshape the data for PCA
pTime('Reshaping Data',t1);
%reshape
testArray2=reshape(testArray,size(testArray,1),size(testArray,2)*size(testArray,3));

%% Flatten patterns properly
%flatten each pattern properly - this is important to get right
testmean=mean(testArray2,1); %calc the mean of each pattern
teststd=std(testArray2,0,1); %calc the std of each pattern
%flatten pattern by pattern
testArray2=(testArray2-repmat(testmean,[size(testArray2,1),1]))./repmat(teststd,[size(testArray2,1),1]);