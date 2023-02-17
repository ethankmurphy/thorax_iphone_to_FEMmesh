%-------------------------------------------------------------------------------
%
%   construct_gmsh_compound_mesh
%
% This code constructs a 3D mesh using gmsh. It relies on an input surface
% triangulation, h-values describing spatially dependent desired tetra
% sizes, and indices for main domain (with subdomain cut out and the
% subdomain. The indices are associated with triangles, which can use the
% same triangles/overlapping triangles. Both sets of triangles have to
% form a water-tight surface. 
% 
% Input Variable:
%          sm_n - Surface mesh nodes 
%          sm_t - Surface mesh triangles
%         hvals - H-values (metric for triangle/tetrahedra size at each
%                 node
%       is_main - Indices of triangles associated with the main domain with
%                 the subdomain cut out.
%         is_sd - Indices of triangles associateed with the subdomain
%   rungmsh_str - path to gmsh
%       
%      Variables:
%          sm_e - Surface mesh edges
%
% Error checking:
%       1. We check that each set of triangle indices (main and subdomain)
%          result in water tight surfaces.
%       2. The union of the two regions is the whole thing (simple check,
%          are all the triangles used)
%       3. (not implemented) The intersection of subregions is empty
% 
%-------------------------------------------------------------------------------
function construct_gmsh_standard_mesh(sm_n,sm_t,hvals,rungmsh_str,geo_str)

if nargin < 5
    geo_str = 'tmp_mesh';%.msh file gmsh makes
end

%-------------------------------------------------------------------------------
% Construct the edges
TR   = triangulation(sm_t, sm_n);%inbuilt - used for convinence
sm_e = edges(TR);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Error checking: Check the surface water tightness
TR = triangulation(sm_t, sm_n);
F  = freeBoundary(TR);
if size(F,1)
    figure
    hold on
    trisurf(sm_t,sm_n(:,1),sm_n(:,2),sm_n(:,3), ...
        'FaceColor','cyan','FaceAlpha', 0.2);
    ks = unique(F(:));
    plot3(sm_n(ks,1),sm_n(ks,2),sm_n(ks,3),'.g','markersize',16)
    error('Surface is not water tight, needs to be fixed')
else
    disp('Water tight, good sign')
end

%-------------------------------------------------------------------------------
% Write the main box nodes to a geo-file
write_nodes2geo('mainbox',sm_n,hvals,0);%makes geo files from matlab

%-------------------------------------------------------------------------------
% Write the all the boundary edges to a geo-file
write_edges2geo('mainbox',sm_e,[0 0]);

%-------------------------------------------------------------------------------
% Write the triangles of the main object to the geo-file
write_tris2geo('mainbox',sm_t,sm_e,[0 0]);

%-------------------------------------------------------------------------------
% Write volume to a geo file
write_vol2geo('mainbox',1:size(sm_t,1),0);

%-------------------------------------------------------------------------------
% Check for the main runner geo file. If its not in the current directory,
% then make it
if filechecker('.','cyl_in_box.geo') == 0
    make_main_geo_runner();
end
%-------------------------------------------------------------------------------
% Check for the temporary .msh and remove it prior to starting the new mesh
if filechecker('gmshes',[geo_str,'.msh']) == 1
    dos(['rm gmshes/',geo_str,'.msh'])
end


%-------------------------------------------------------------------------------
% Construct the gmesh mesh: Run and save a mesh in gmsh
run_str = [                           ...
    rungmsh_str                     , ...
    'cyl_in_box.geo'                , ...
    ' -3 '                          , ...%make 3d volume    
    ' -optimize_threshold 0.6 -format msh22 '     , ...%optimise
    '-o gmshes/',geo_str            , ...%outputfilename
    '.msh'];%command line for gmsh

dos(run_str)%dos sends to commandline oustide matlab


%-------------------------------------------------------------------------------
% 
% Make the main geo file
%
%-------------------------------------------------------------------------------
function make_main_geo_runner()

%-------------------------------------------------------------------------------
% Open the file
fid = fopen('cyl_in_box.geo','w');

%-------------------------------------------------------------------------------
fprintf(fid,'/****************************************/\n');
fprintf(fid,'// Construct the nodes of the boundary and the interior probe/cylinder\n');
fprintf(fid,'Include "mainbox_nodes.geo";\n\n');

fprintf(fid,'/****************************************/\n');
fprintf(fid,'// Construct the edges of the boundary and the interior probe/cylinder\n');
fprintf(fid,'Include "mainbox_edges.geo";\n\n');

fprintf(fid,'/****************************************/\n');
fprintf(fid,'// Construct the triangles or patches\n');
fprintf(fid,'Include "mainbox_tris.geo";\n');

fprintf(fid,'/****************************************/\n');
fprintf(fid,'// Add in tetrahedra representing the coarse elements, which also combines all\n');
fprintf(fid,'// the elements\n');
fprintf(fid,'Include "mainbox_vol.geo";\n\n');

fprintf(fid,'check = 0;\n\n');

%-------------------------------------------------------------------------------
% Close the file
fclose(fid);