%--------------------------------------------------------------------------
%
% Unwrap the boundary and make a 2D mesh with encoded electrodes, then 
% wrap it back into 3D. The general steps are as follows:
%   1. Rotate the mesh so that the +x-axis is aligned with the largest gap
%      between electrodes. This helps ensure there is no issues with
%      unwrapping and wrapping. 
%   2. Convert all the points to a 2D frame of arc length vs height
%   3. Make the mesh with encoded electrodes using distmesh. First we need
%      to set boundary points, then we just call distmesh
%   4. Check that there are matching sets of points on the left and right
%      side of the domain so, we can actually stitch it back together
%   5. Stitch things back together and bring things to 3D
%   6. Do some adjustments to the electrodes, so that the area is nearly
%      what we want them to be while also being encoding/comforming to the
%      surface
%   7. Unrotate the meshes to their original state and end
%
%--------------------------------------------------------------------------
function [t,bnd_f,el3D_bndpts,area_errs] = unwrap_bnd_v3(el3D_bndpts,delec,bmsh,elecs,bndfit,elns,dbg_flg)


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   1. Rotate the mesh so that the +x-axis is aligned with the largest gap
%      between electrodes. This helps ensure there is no issues with
%      unwrapping and wrapping. 
%--------------------------------------------------------------------------
% Calculate the minimum distance between electrode boundary points
minds  =  1000*ones(length(el3D_bndpts),3);
minzel =  100000;
maxzel = -100000;
for n = 1:length(el3D_bndpts)
    if n == length(el3D_bndpts)
        k = 1;
    else
        k = n+1;
    end
    minds(n,:) = [n k min(eval_closest_pdist(el3D_bndpts{n},el3D_bndpts{k},0,1))];
    minzel     = min([minzel; el3D_bndpts{n}(:,3)]);
    maxzel     = max([maxzel; el3D_bndpts{n}(:,3)]);
end
[tmp,i] = max(minds(:,3));
telec   = atan2(elecs(:,2),elecs(:,1));
tref    = mean(telec(minds(i,1:2)));
if abs(telec(minds(i,1)) - telec(minds(i,2)) ) > pi/4
    error('stop')
end
%--------------------------------------------------------------------------
% Rotate the whole mesh so the center of this largest space is right on
% x-axis
R0 = [cos(tref) sin(tref) 0; -sin(tref) cos(tref) 0; 0 0 1];
bmsh.node = (R0*(bmsh.node'))';
elecs     = (R0*(elecs'))';
for n = 1:length(el3D_bndpts)
    el3D_bndpts{n} = (R0*(el3D_bndpts{n}'))';
end
if dbg_flg == 1
    figure
    hold on
    trisurf(bmsh.tri,bmsh.node(:,1),bmsh.node(:,2),bmsh.node(:,3), 'Facecolor','cyan','FaceAlpha',0.8,'linestyle','none');
    axis equal;
    for n = 1:length(el3D_bndpts)
        elps = el3D_bndpts{n};
        plot3(elps(:,1),elps(:,2),elps(:,3),'ok','markersize',2)
        text(1.1*mean(elps(:,1)),1.1*mean(elps(:,2)),mean(elps(:,3)),num2str(n))
    end
    lbl_fmt_fig('X','Y','Rotate so the maximum spacing is on the +x-axis',[],'Z',16)
    view([10 6])
    minds(i,:)
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   2. Convert all the points to a 2D frame of arc length vs height
%--------------------------------------------------------------------------
minz = min(bmsh.node(:,3));
maxz = max(bmsh.node(:,3));
%--------------------------------------------------------------------------
% Perform a smooth fit to the rotated surface
[pfits, rbf4ps]  = four_cylfit(bmsh.node,bndfit.Nf,bndfit.Nzc,bndfit.sigz*10);
% figure;hold on
% trisurf(bmsh.tri,bmsh.node(:,1),bmsh.node(:,2),bmsh.node(:,3),'facecolor','cyan','linestyle','none','facealpha',0.3)
% plot3(pfits(:,1),pfits(:,2),pfits(:,3),'.k','markersize',12)

%--------------------------------------------------------------------------
% Convert the points to cylindrical coordinates
ts       = linspace(0,2*pi,64)';
zs       = mean(elecs(:,3))+(0*ts);
bfitxyzs = eval_four_cylfit([ts zs],rbf4ps);
rxys     = sqrt( bfitxyzs(:,1).^2 + bfitxyzs(:,2).^2 );
%--------------------------------------------------------------------------
% Convert the angles to a distance by assuming a cylinder with the mean
% radius of all the radii
[s,transobj] = fourb_trans(rxys,ts,[],1);
%--------------------------------------------------------------------------
% Transform the 3D electrode points points to the 2D points defined as the
% arc length and the height
tels  = [];
telcs = [];
zels  = [];
sels  = [];
for n = 1:length(el3D_bndpts)
    tmps        = atan2(el3D_bndpts{n}(:,2),el3D_bndpts{n}(:,1));
    isneg       = find(tmps < 0);
    tmps(isneg) = tmps(isneg) + 2*pi;
    %----------------------------------------------------------------------
    tels       = [tels;  tmps];
    telcs      = [telcs; mean(tmps)];
    zels       = [zels;  el3D_bndpts{n}(:,3)];
    seltmp     = fourb_trans([],tmps,transobj,2);
    sels       = [sels; seltmp];
    elpts2D{n} = [seltmp el3D_bndpts{n}(:,3)];
end
elpts = [sels zels];


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   3. Make the mesh with encoded electrodes using distmesh. First we need
%      to set boundary points, then we just call distmesh
%--------------------------------------------------------------------------
% Add the boundary points: Add enough points on the left and right side so
% that these represent all the end points. This will allow us to easily 
% stitch the domain back together
tstt  = 0;
tend  = tstt     + 2*pi;
%-----------------------------
sstt  = 0;
send  = fourb_trans([], tend,transobj,2);
%-----------------------------
% Parameters here are reasonable for adult thorax domains assuming units of
% mm. If there is a big change, then these values need to be changed. 
% distance/size = number, or number*size = distance
zspc  = 1;
ntop  = round(abs(   maxz - (maxzel+zspc) )/4);
nbot  = round(abs(   minz - (minzel-zspc) )/4);
nmid  = round(abs(maxzel - minzel)/2);
nsstp = round(abs(send)/10); % was 10
zs    = unique([ ...
    linspace(       minz,minzel-zspc,nbot)'; ...    % coarse bottom
    linspace(minzel-zspc,maxzel+zspc,nmid)'; ...    % fine middle
    linspace(maxzel+zspc,maxz,ntop)']);             % coarse top      
nz    = length(zs);   
bbnd  = [ ...
    linspace(   0,send,nsstp)' minz*ones(nsstp,1);                  % Bottom side
    linspace(   0,send,nsstp)' maxz*ones(nsstp,1);                  % Top side
    sstt*ones(nz,1) zs
    send*ones(nz,1) zs];
bbnd = unique(bbnd,'rows');
disp(['distmesh input: zspc=',num2str(zspc),', N_horz=',num2str(nsstp),', N_bndps=',num2str(size(bbnd,1))])
%--------------------------------------------------------------------------
% Make sure the boundaries are arranged in an CCW fashion
ptmp     = bbnd - repmat(mean(bbnd,1),size(bbnd,1),1);
ts       = atan2(ptmp(:,2),ptmp(:,1));
is       = find(ts<pi/2);
ts(is)   = ts(is)+2*pi;
[tmp,is] = sort(ts);
bbnd     = bbnd(is,:);
%--------------------------------------------------------------------------
% Construct the mesh
[p,t,hvals] = construct_distmesh_polygon_innerfixed(bbnd,elpts,dbg_flg,1.2);
if dbg_flg == 1
    hold on
    plot(elpts(:,1),elpts(:,2),'xg')
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   4. Check that there are matching sets of points on the left and right
%      side of the domain so, we can actually stitch it back together
%--------------------------------------------------------------------------
% Verify with a plot
if dbg_flg == 1
    figure
    subplot(1,2,1)
    hold on
    is  = find( abs(p(:,1) - min(p(:,1)))<0.001 );
    plot(p(is,2),'ok')
    is  = find(bbnd(:,1) == min(bbnd(:,1)) );
    plot(bbnd(is,2),'xr')
    legend('Distmesh p','Input fixed p')
    title('Left side')
    subplot(1,2,2)
    hold on
    is  = find(abs(p(:,1) - max(p(:,1)))<0.001 );
    plot(p(is,2),'ok')
    is  = find(bbnd(:,1) == max(bbnd(:,1)) );
    plot(bbnd(is,2),'xr')
    title('Right side')
end
%--------------------------------------------------------------------------
% Connect the end points: For all the far right points redefine these as
% the analogous far left points
il  = find( abs(p(:,1) - min(p(:,1)))<0.001 );
ir  = find( abs(p(:,1) - max(p(:,1)))<0.001 );
ir  = ir(end:-1:1);
% save testdata
if abs(length(il) - length(ir)) > 0
    % Trying again with a larger h-increase factor
    [p,t,hvals] = construct_distmesh_polygon_innerfixed(bbnd,elpts,dbg_flg,1.4);
    il  = find( abs(p(:,1) - min(p(:,1)))<0.001 );
    ir  = find( abs(p(:,1) - max(p(:,1)))<0.001 );
    ir  = ir(end:-1:1);
end
if norm( p(il,2) - p(ir,2) ) > 0.001
    ir  = ir(end:-1:1); % Try switching it again!!
    if norm( p(il,2) - p(ir,2) ) > 0.001
        norm( p(il,2) - p(ir,2) )
        error('points don''t match up')
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   5. Stitch things back together and bring things to 3D
%--------------------------------------------------------------------------
for n = 1:length(ir)
    for k = 1:3
        i0      = find( t(:,k) == ir(n));
        t(i0,k) = il(n);
    end
end
%--------------------------------------------------------------------------
% Remove unused nodes
[t,p] = remove_unused_nodes(t,p);
%--------------------------------------------------------------------------
% Bring back to a chest-shaped boundary
ts_f  = fourb_trans([],p(:,1),transobj,-1);
bnd_f = eval_four_cylfit([ts_f p(:,2)],rbf4ps);
%--------------------------------------------------------------------------
for n = 1:length(el3D_bndpts)
    elpts = elpts2D{n};
    %----------------------------------------------------------------------
    % Bring back to a chest-shaped boundary
    ts_el      = fourb_trans([],elpts(:,1),transobj,-1);
    elpts3D{n} = eval_four_cylfit([ts_el elpts(:,2)],rbf4ps);
end
%---------------------------------------------------------------------------
% Remove triangles that connect the ends
rt = [];
for n = 1:size(t)
    d1 = norm( p(t(n,1),:) - p(t(n,2),:));
    d2 = norm( p(t(n,2),:) - p(t(n,3),:));
    d3 = norm( p(t(n,1),:) - p(t(n,3),:));
    if max([d1 d2 d3]) > 10
        rt = [rt; n];
    end
end
tcut = t;
tcut(rt,:) = [];



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   6. Do some adjustments to the electrodes, so that the area is nearly
%      what we want them to be while also being encoding/comforming to the
%      surface
%--------------------------------------------------------------------------
area_errs = zeros(length(elpts2D),2);
for n = 1:length(elpts2D)
    %----------------------------------------------------------------------
    % Adjust the electrode points until the area is equal to our
    % requisted electrode area
    bestfac = fminbnd(@(fac) elec_pts_fit_area(fac,n,delec,p,elecs,elpts2D,elns,bnd_f,tcut),0.60,1.40);
    % bestfac
    %----------------------------------------------------------------------
    % Get the boundary points
    e3Dis = zeros(length(el3D_bndpts{n}),1);
    for k = 1:length(el3D_bndpts{n})
        [tmp,e3Dis(k)] = min( (bnd_f(:,1)-el3D_bndpts{n}(k,1)).^2 + (bnd_f(:,2)-el3D_bndpts{n}(k,2)).^2 + (bnd_f(:,3)-el3D_bndpts{n}(k,3)).^2);
    end
    
    %----------------------------------------------------------------------
    % Keep the best nodes
    [area_err,aTRI_E,bnd_f] = elec_pts_fit_area(bestfac,n,delec,p,elecs,elpts2D,elns,bnd_f,tcut);
    if dbg_flg == 1
        disp(['E',num2str(n),': Area error = ',num2str(area_err),', fac = ',num2str(bestfac)])
    end
    area_errs(n,1) = area_err;
    area_errs(n,2) = bestfac;
    %----------------------------------------------------------------------
    % Update the 2D points
    tmps       = atan2(bnd_f(e3Dis,2),bnd_f(e3Dis,1));
    %----------------------------------------------------------------------
    % If any angles are wrapped around fix this
    D = comp_pairwise_distmat(tmps);
    if max(D(:)) > pi/2
        if mean(tmps) > pi/2
            is = find(tmps < pi/2);
            tmps(is) = tmps(is) + 2*pi;
        else
            is = find(tmps > pi/2);
            tmps(is) = tmps(is) - 2*pi;
        end
    end
    seltmp     = fourb_trans([],tmps,transobj,2);
    elpts2D{n} = [seltmp bnd_f(e3Dis,3)];

    %----------------------------------------------------------------------
    % Update the 3D points to include only points on the mesh
    el3D_bndpts{n} = bnd_f(e3Dis,:);

end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   7. Unrotate the meshes to their original state and end
bnd_f     = ((R0')*(bnd_f'))';
elecs     = ((R0')*(elecs'))';
for n = 1:length(el3D_bndpts)
    el3D_bndpts{n} = ((R0')*(el3D_bndpts{n}'))';
end


%--------------------------------------------------------------------------
% Plot the chest domain surface and overlay the electrode points
if dbg_flg == 1
    figure
    hold on
    trisurf(t,bnd_f(:,1),bnd_f(:,2),bnd_f(:,3), 'Facecolor','cyan','FaceAlpha',0.8); axis equal;
    for n = 1:length(el3D_bndpts)
        elpts = elpts2D{n};
        %------------------------------------------------------------------
        % Bring back to a chest-shaped boundary
        ts_el   = fourb_trans([],elpts(:,1),transobj,-1);
        bnd_el  = eval_four_cylfit([ts_el elpts(:,2)],rbf4ps);
        norm( bnd_el - el3D_bndpts{n})
        %------------------------------------------------------------------
        plot3(bnd_el(:,1),bnd_el(:,2),bnd_el(:,3),'-or')
    end
    lbl_fmt_fig('X','Y','Wrapped Distmesh surface',[],'Z',16)
    view([20 30])
    saveas(gcf,'step3_wrapped_3D_distmsh','png')
end


