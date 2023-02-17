%-------------------------------------------------------------------------------
%
% Adjust the electrode points until the area is equal to our
% requisted electrode area
%
%-------------------------------------------------------------------------------
function [area_err,aTRI_E,bnd_f_tmp] = elec_pts_fit_area(fac,n,delec,pts,elc3D,elpts2D,elns,bnd_f,t)

%-------------------------------------------------------------------------------
% Adjust the nth electrode points by a factor (in the radial
% direction
% * Find the nodes of the nth electrode (including the boundary
%   nodes
% * Find 2 orthonormal vectors that are perpendicular to the electrode
%   normal vector
% * Adjust the points in the radial direction
%-------------------------------------------------------------------------------
[IN,ON] = inpolygon(pts(:,1),pts(:,2),elpts2D{n}(:,1),elpts2D{n}(:,2));
IN      = find( (IN == 1) | (ON == 1));
%-------------------------------------------------------------------------------
% Get two orthonormal vectors in the plane
nvec      = elns(n,:);
[uv1,uv2] = get_two_orthvecs(nvec);
%-------------------------------------------------------------------------------
% Stretch (or shrink) points in the radial direction
bnd_f_tmp = bnd_f;
for k = 1:length(IN)
    rel3Dpt            = bnd_f_tmp(IN(k),:) - elc3D(n,:);
    relprj3D           = fac*dot(rel3Dpt,uv1)*uv1 + fac*dot(rel3Dpt,uv2)*uv2 + dot(rel3Dpt,nvec)*nvec;
    prj3Dpt            = relprj3D       + elc3D(n,:);
    bnd_f_tmp(IN(k),:) = prj3Dpt;
end

%-------------------------------------------------------------------------------
% Calculate the area of the nth electrode
% * Get the triangles involved in the nth electrode
% * Calculate all the areas of the triangles in 3D
% * Sum the areas of the nth electrode
%-------------------------------------------------------------------------------
pts2Dt  = 1/3*(pts(t(:,1),:)+pts(t(:,2),:)+pts(t(:,3),:));
INtri   = find( inpolygon(pts2Dt(:,1),pts2Dt(:,2),elpts2D{n}(:,1),elpts2D{n}(:,2)) == 1);
areaTRI = calc_TRI_area(t,bnd_f_tmp);
aTRI_E  = sum( areaTRI(INtri) );

%-------------------------------------------------------------------------------
area_err = abs( pi*(delec/2)^2 - aTRI_E);