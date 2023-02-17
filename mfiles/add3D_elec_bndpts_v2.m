%--------------------------------------------------------------------------
%
% Construct 3D electrodes in the plane of the normal vector
% 
%--------------------------------------------------------------------------
function [el3D_bndpts,elns] = add3D_elec_bndpts_v2(elc3D,delec,bmsh,dbg_flg)

%--------------------------------------------------------------------------
% Calculate all the normal vectors and make sure they are pointing towards
% the surface
nvcs   = calc_nvecs(bmsh.node,bmsh.tri);
tcnts  = get_tcs(bmsh.node,bmsh.tri);
isflip = find( sum(nvcs.*tcnts,2) < 0 );
nvcs(isflip,:) = -nvcs(isflip,:);

%--------------------------------------------------------------------------
% Calculate the normal vector at each electrode by using a weighted 
% average of nearby triangles
elns  = zeros(size(elc3D,1),3);
for n = 1:size(elc3D,1)
    %----------------------------------------------------------------------
    % Calculate distances from the electrode
    dists = sqrt((tcnts(:,1)-elc3D(n,1)).^2 + (tcnts(:,2)-elc3D(n,2)).^2 + (tcnts(:,3)-elc3D(n,3)).^2);
    
    %----------------------------------------------------------------------
    % Caclulate a normal vector at the electrode point. Use a weighted 
    % average of the centers that are within a diameter of the elctrode 
    % center. Unless there are no centers, then using the closest
    is = find( dists < delec/4 );
    if isempty(is) == 1
        [tmp,i0] = min(dists);
        ntmp     = nvcs(i0,:);
    else
        ntmp     = sum(repmat(dists(is),1,3).*nvcs(is,:),1)./repmat(sum(dists(is)),1,3);
    end
    elns(n,:) = ntmp;
end

%--------------------------------------------------------------------------
% Sample circular boundary points centered at each electrode with distance
% delec from it and normal to the electrode direction
Np       = 16;
ts       = linspace(0,2*pi,Np+1)';
ts       = ts(1:end-1);
for n = 1:size(elc3D,1)
    %----------------------------------------------------------------------
    % Get two perpendicular vectors that are perpendicular to the normal
    % vector
    %----------------------------------------------------------------------
    % Unit vector 1
    vv1 = rand(1,3);
    nv  = elns(n,:);
    uv1 = cross(nv,vv1);
    uv1 = uv1 / norm(uv1);
    %----------------------------------------------------------------------
    % Unit vector 2
    uv2 = cross(nv,uv1);
    uv2 = uv2 / norm(uv2);
    
    %----------------------------------------------------------------------
    elps = repmat(elc3D(n,:),Np,1) +   ...
        repmat(delec/2*cos(ts),1,3).*repmat(uv1,Np,1) + ...
        repmat(delec/2*sin(ts),1,3).*repmat(uv2,Np,1);
    
    %----------------------------------------------------------------------
    % Record the electrode points
    el3D_bndpts{n} = elps;
        
end

if dbg_flg == 1
    [(1:size(elc3D,1))' round(atan2(elc3D(:,2),elc3D(:,1))*180/pi,1)]
    
    figure
    hold on
    trisurf(bmsh.tri,bmsh.node(:,1),bmsh.node(:,2),bmsh.node(:,3), 'Facecolor','cyan','FaceAlpha',0.8,'linestyle','none');
    axis equal;
    for n = 1:size(elc3D,1)
        elps = el3D_bndpts{n};
        plot3(elps(:,1),elps(:,2),elps(:,3),'ok','markersize',2)
        text(1.1*mean(elps(:,1)),1.1*mean(elps(:,2)),mean(elps(:,3)),num2str(n))
    end    
    
    lbl_fmt_fig('X','Y','Electrode Boundaries',[],'Z',16)
    view([10 6])
    % saveas(gcf,'figs/step2_add_electrode_bounds','png')
end
