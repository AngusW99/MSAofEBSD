%% Read the EDX spectra
pattern_number=1;
[EDSData_cor,EDSData_raw ] = bReadEDX(EBSPData,pattern_number,channum);

EDSData_cor_map=zeros(Data_InputMap.ypts,Data_InputMap.xpts,channum);
EDSData_raw_map=zeros(Data_InputMap.ypts,Data_InputMap.xpts,channum);

Mean_cor=zeros(Data_InputMap.ypts,Data_InputMap.xpts);
Mean_raw=zeros(Data_InputMap.ypts,Data_InputMap.xpts);

pTime('Loading EDX Data',t1);
for xi=1:Data_InputMap.xpts
    parfor yi=1:Data_InputMap.ypts
        [EDSData_cor_map(yi,xi,:),EDSData_raw_map(yi,xi,:) ] = bReadEDX(EBSPData,Data_InputMap.PMap(yi,xi),channum);
        Mean_cor(yi,xi)=mean(EDSData_cor_map(yi,xi,:));
        Mean_raw(yi,xi)=mean(EDSData_raw_map(yi,xi,:));
    end
end

%reshape the array
EDSData_cor_map2=permute(EDSData_cor_map,[3,1,2]);
EDSData_cor_map2=reshape(EDSData_cor_map2(1:2048/2,:,:),[2048/2,size(EDSData_cor_map2,2)*size(EDSData_cor_map2,3)]);