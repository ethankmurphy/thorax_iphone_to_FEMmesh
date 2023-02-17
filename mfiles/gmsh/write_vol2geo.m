%-------------------------------------------------------------------------------
%
% Write the boundary points to a geo file
%
%------------------------------------------------------------------------------- 
function write_vol2geo(prfx,vol_inds,ios)   

%-------------------------------------------------------------------------------
% Open the geo file
fid = fopen([prfx,'_vol.geo'],'w');

%-------------------------------------------------------------------------------
% Loop through the points    
fprintf(fid,'Surface Loop(%i) = {',1+ios);
for n = 1:length(vol_inds)-1
    fprintf(fid,'%i,',vol_inds(n)+ios);
end
fprintf(fid,'%i};\n',vol_inds(end)+ios);
fprintf(fid, ...
    'Volume(%i) = {%i};\n'              , ...
    1+ios                               , ...
    1+ios );

fclose(fid);