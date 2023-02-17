%-------------------------------------------------------------------------------
%
% Write the boundary points to a geo file
%
%------------------------------------------------------------------------------- 
function write_tris2geo(prfx,tris,Es,ios)   

%-------------------------------------------------------------------------------
% Open the geo file
fid = fopen([prfx,'_tris.geo'],'w');
fprintf(fid,'p    = %i;\n',ios(2));
%-------------------------------------------------------------------------------
% Loop through the points
for n = 1:size(tris,1)
    %---------------------------------------------------------------------------
    ll = zeros(1,3);
    for k = 1:3
        if k < 3         
            j = find( (Es(:,1) == tris(n,k)) & (Es(:,2) == tris(n,k+1)) )+ios(1);
        else
            j = find( (Es(:,1) == tris(n,k)) & (Es(:,2) == tris(n,  1)) )+ios(1);
        end
        if isempty(j) == 1
            if k < 3
                j = -(find( (Es(:,2) == tris(n,k)) & (Es(:,1) == tris(n,k+1)) )+ios(1));
            else
                j = -(find( (Es(:,2) == tris(n,k)) & (Es(:,1) == tris(n,1)) )+ios(1));
            end
        end
        ll(k) = j;        
    end
    %---------------------------------------------------------------------------
    fprintf(fid,'p = p+1;\n');
    fprintf(fid, ...
        'Line Loop(p) = {%i,%i,%i};\n', ...
        ll(1), ll(2), ll(3) );    
    fprintf(fid,'Plane Surface(p) = {p};\n');    
    % Line Loop(llid) = {  ios+ls_1,  ios+ls_2,  -(ios+ls_3),  -(ios+ls_4)};   
    % sid = sid+1;//news;
    % //Plane Surface(sid) = {llid};
end
% if ios(2) == 0
%     fprintf(fid,'Surface Loop( 1 ) = { 1:p};\n');
%     fprintf(fid,'Volume(1) = {1};\n');
%     % else
%     %     fprintf(fid,'Surface {(%i +1):p } In Volume {1};\n',ios);
% end

%-------------------------------------------------------------------------------
% Close the file
fclose(fid);