%% Read in the data
pTime('Reading in h5 data',t1);
% open HDF5 file and read in some data on EBSD map
% InputUser.HDF5_folder=BaseFolder;
% InputUser.HDF5_file=HDF5_filename;
[ MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
%read the map & convert to area data
[Data_InputMap] = EBSD_Map(MapData,MicroscopeData);
