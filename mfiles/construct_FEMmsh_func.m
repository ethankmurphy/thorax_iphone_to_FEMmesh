%--------------------------------------------------------------------------
% 
%
%
%--------------------------------------------------------------------------
function [msh,eareas] = construct_FEMmsh_func(bp,t,el3D_bndpts,delec,hbnd,hel,elecs,rungmsh_str,dbg_flg)


%--------------------------------------------------------------------------
% Set the h-values
hvals      = hbnd*ones(size(bp,1),1);
iel        = [];
for n1 = 1:length(el3D_bndpts)
    for n2 = 1:size(el3D_bndpts{n1},1)
        iel = [iel; find( sqrt(sum( (bp - el3D_bndpts{n1}(n2,:)).^2,2)) < delec )];
    end
end
hvals(iel) = hel;
if dbg_flg == 1
    figure; hold on
    trisurf(t,bp(:,1),bp(:,2),bp(:,3),'facecolor','cyan', ... 'linestyle','none', ...
        'facealpha',0.5)
    plot3(bp(iel,1),bp(iel,2),bp(iel,3),'.m','markersize',16)
    axis equal
    view([10 6])
end

%--------------------------------------------------------------------------
% Construct the mesh
construct_gmsh_standard_mesh(bp,t,hvals,rungmsh_str);

%--------------------------------------------------------------------------
% Extract the nodes and elements from the gmsh mesh file
msh = read_gmsh_mesh_v3('gmshes/tmp_mesh',0);

%--------------------------------------------------------------------------
% 2.b.  Determine the face elements using NDRM's face mex face
%     finding algorithm
disp('Finding faces')
TR       = triangulation(msh.elem, msh.node);
msh.face = freeBoundary(TR);

%--------------------------------------------------------------------------
% 3. Find the electrodes and their respective faces
% Make sure elec is a cell array of size (number of elecs) x 1, and
% each elec component (elec{n}) is size 1 x (Number of face elements)
fcnts = 1/3*(msh.node(msh.face(:,1),:) + msh.node(msh.face(:,2),:) + msh.node(msh.face(:,3),:));

%------------------------------------------------------------------
for n = 1:length(el3D_bndpts)
    % Get a best fitting plane to the current electrode       
    bps = el3D_bndpts{n};
    fcs_tmp = fcnts;
    %-----------------------
    [nvec,cent,long_ax,shrt_ax] = get_nrmal_vec(bps);
    % Transform points to 2D
    bps_2D = [ ...
        sum( (bps-repmat(cent,size(bps,1),1)).*repmat(long_ax',size(bps,1),1),2) ...
        sum( (bps-repmat(cent,size(bps,1),1)).*repmat(shrt_ax',size(bps,1),1),2) ...
        sum( (bps-repmat(cent,size(bps,1),1)).*repmat(nvec',size(bps,1),1),2)];
    % face centers
    fcs_2D = [ ...
        sum( (fcs_tmp-repmat(cent,size(fcnts,1),1)).*repmat(long_ax',size(fcnts,1),1),2) ...
        sum( (fcs_tmp-repmat(cent,size(fcnts,1),1)).*repmat(shrt_ax',size(fcnts,1),1),2) ...
        sum( (fcs_tmp-repmat(cent,size(fcnts,1),1)).*repmat(nvec',size(fcnts,1),1),2) ];
    % Make sure the electrode points form a polygon
    av_bps_2D = mean(bps_2D,1);
    thts      = atan2(bps_2D(:,2)-av_bps_2D(2),bps_2D(:,1)-av_bps_2D(1));
    [tmp,sid] = sort(thts);
    bps_2D    = bps_2D(sid,:);
    % Find points in the polygon
    is = find( (inpolygon(fcs_2D(:,1),fcs_2D(:,2),1.02*bps_2D(:,1),1.02*bps_2D(:,2)) == 1) & ...
        (abs(fcs_2D(:,3)) < 10) );
    %-----------------------
    elec{n,1} = is';
end
msh.elec = elec;

%--------------------------------------------------------------------------
% Calculate the electrode areas
eareas = zeros(length(elec),1);
for n = 1:length(elec)
    eareas(n) = sum( calc_TRI_area(msh.face(elec{n},:),msh.node,0) );
end


%--------------------------------------------------------------------------
% Convert to meters
msh.node = msh.node/1000;
