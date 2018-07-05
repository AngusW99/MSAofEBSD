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