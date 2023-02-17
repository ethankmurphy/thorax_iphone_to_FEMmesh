%-------------------------------------------------------------------------------
%
% Construct a 2D mesh over an inner and outer polygon domain. The h-values
% of the mesh should match the edge spacing of the inner and outer
% boundaries. 
%
%-------------------------------------------------------------------------------
function [p,t] = construct_distmesh_polygon_nofixed(pout,dbg_flg,inchval)

%-------------------------------------------------------------------------------
% Make sure the inner and outer boundaries are arranged in an CCW fashion
% Outer:
ptmp     = pout - repmat(mean(pout,1),size(pout,1),1);
ts       = atan2(ptmp(:,2),ptmp(:,1));
[tmp,is] = sort(ts);
pout     = pout(is,:);

% figure
% plot(pout(:,1),pout(:,2),'-xg')
%-------------------------------------------------------------------------------
% H-values: Calculate the spacing of points on the inside portion and the 
% outside portion
% Outer:
D      = comp_pairwise_distmat(pout);       
minDs  = sort(D,2,'ascend');
if nargin < 3
    hout   = 1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
else
    hout   = inchval*1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
end

%-------------------------------------------------------------------------------
% Distmesh
%-------------------------------------------------------------------------------
% Construct a distance function with the drill bit cut out. In this
% function the drillbit is much more complicated than a circle
fd=@(p) dpoly(p,[pout; pout(1,:)]);
%-------------------------------------------------------------------------------
% Set the bounding box and the fixed points
fpts = pout;
bbx  = [min(pout,[],1); max(pout,[],1)];  
%-------------------------------------------------------------------------------
figure
[p,t]=distmesh2d(fd,@huniform,hout,bbx,fpts);%not relates to p above. p is n x 3, t is triangulation of dimen m x 3

%-------------------------------------------------------------------------------
if dbg_flg == 1 
    hold on
    plot(fpts(:,1),fpts(:,2),'dg','linewidth',2,'markersize',3)
    lbl_fmt_fig('X','Y','')
    axis equal
    
end

