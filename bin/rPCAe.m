%% do PCA on the EDS

[coeffe,scoree,latente,~,explainede]=pca(EDSData_cor_map2, 'Centered',false, 'NumComponents',NumPCAe);

%% rotate factors
[coeffVMe, RotVMe] = rotatefactors(coeffe(:,1:NumPCAe),'Method','varimax', 'Maxit', 5000, 'Reltol', 0.5*sqrt(eps));
scoreVMe=scoree(:,1:NumPCAe)*RotVMe;
%coeff = maps
%scores = patterns