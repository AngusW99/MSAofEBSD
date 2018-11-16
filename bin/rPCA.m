%%  run PCA analysis to split data into basis patterns (eigenvectors) and signal strength in map (loadings)
%   this uses a user input choice on number of factors to retain (NumPCA)

%coeff = maps
%scores = patterns

%% perform the PCA
pTime('PCA normal',t1);
[coeff,score,latent,~,explained]=pca(testArray2, 'Centered',false,'NumComponents',NumPCA);

%% Rotate the factors to improve contrast in the factors (i.e. make the indexable)

pTime('Rotating Factors',t1);
%rotate the factors
[coeffVM, RotVM] = rotatefactors(coeff(:,1:NumPCA),'Normalize','on','Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
scoreVM=score(:,1:NumPCA)*RotVM;

%% Sort the factor loadings to make them positive

% sort the factor map
PCA_VM_Pat_n=zeros(PatSizeW,PatSizeH,NumPCA); %pattern
PCA_VM_Map_n=zeros([Data_InputMap.ypts Data_InputMap.xpts NumPCA]); %load maps
PCA_VM_num=zeros([Data_InputMap.ypts Data_InputMap.xpts]); 
PCA_VM_val=zeros([Data_InputMap.ypts Data_InputMap.xpts]); %pattern to use

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