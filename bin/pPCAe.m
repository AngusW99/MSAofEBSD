%plot the EDS based PCAs
figure;

for n=1:NumPCAe
    subplot(PCAGride(1),PCAGride(2),n);
    iqr_1=iqr(coeffVMe(:,n));
    med_1=median(coeffVMe(:,n));
    crange=med_1+[-iqr_1/2 iqr_1/2];
    imagesc(reshape(coeffVMe(:,n),[size(EDSData_cor_map,1),size(EDSData_cor_map,2)])); axis image; caxis(crange);
end

for n=1:NumPCAe
    subplot(PCAGride(1),PCAGride(2),n+3);
    plot(scoree(:,n));
    xlim([0 2048/2]);
end

