function [xi,yi,RefPat_cor,RefPat]=grabEBSP(f1)
%grab an EBSP from a map to inspect
figure(f1);
[xi,yi]=ginput(1);
xi=round(xi); yi=round(yi);
hold on; scatter(xi,yi);

pattern_number=Data_InputMap.PMap(yi,xi);
[ RefPat ] = bReadEBSP(EBSPData,pattern_number);
[ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );

figure; imagesc(RefPat_cor);axis image; axis ij; axis tight; axis off;title(['X = ' int2str(n) ' Y = ' int2str(yi) ])