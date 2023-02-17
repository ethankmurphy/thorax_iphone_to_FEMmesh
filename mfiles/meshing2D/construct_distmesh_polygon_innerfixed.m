%-------------------------------------------------------------------------------
%
% Construct a 2D mesh over an inner and outer polygon domain. The h-values
% of the mesh should match the edge spacing of the inner and outer
% boundaries. 
%
%-------------------------------------------------------------------------------
function [p,t,hvals] = construct_distmesh_polygon_innerfixed(pout,fps,dbg_flg,inchval)

if nargin < 4
    inchval = 1;
elseif nargin < 3
    inchval = 1;
    dbg_flg = 0;
end
%-------------------------------------------------------------------------------
% Make sure the inner and outer boundaries are arranged in an CCW fashion
% Outer:
ptmp     = pout - repmat(mean(pout,1),size(pout,1),1);
ts       = atan2(ptmp(:,2),ptmp(:,1));
is = find(ts<pi/2);
ts(is)   = ts(is)+2*pi;
[tmp,is] = sort(ts);
pout     = pout(is,:);

%-------------------------------------------------------------------------------
% H-values: Calculate the spacing of points on the inside portion and the 
% outside portion
% Inner:
D      = comp_pairwise_distmat(fps);
minDs  = sort(D,2,'ascend');
hin    = inchval*1.2 * sqrt(max( minDs(:,2))) * sin(60*pi/180);
% Outer:
D      = comp_pairwise_distmat(pout);
minDs  = sort(D,2,'ascend');
hout   = inchval*1.1 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
houtv  = inchval*1.1 * sqrt(minDs(:,2)) * sin(60*pi/180);
% Sort the two h-values (small to big)
hmin   = min([hin hout]);
hvals  = [hin hout];

%-------------------------------------------------------------------------------
% Distmesh: Make the 2D mesh with whole
%-------------------------------------------------------------------------------
% Construct a distance function with the drill bit cut out. In this
% function the drillbit is much more complicated than a circle
fd=@(p,binf) dpoly(p,[pout; pout(1,:)]);
%-------------------------------------------------------------------------------
% Set the bounding box and the fixed points
fpts = [fps; pout];
bbx  = [min(pout,[],1); max(pout,[],1)];    
%-------------------------------------------------------------------------------
% Get a convex hull about the fpts, to help provide a max stretch for those
% points that may be interior
K = convhull(fps(:,1),fps(:,2));
% whos
% figure
% hold on
% plot(fps(:,1),fps(:,2),'.g','markersize',12)
% plot(fps(K,1),fps(K,2),'-k')
%-------------------------------------------------------------------------------
% Construct the h-scaling function:
%   The scaling function is basic. For each input point we find the closest
%   Inner point and closest outer point. Then we construct a linear
%   function that changes from hin to hout based on where it is between the
%   two points. 
binf.hin   = hin;
binf.hout  = hout;
binf.houtv = houtv;
binf.pin   = fps;
binf.pout  = pout;
binf.K     = K;
if dbg_flg == 2
    % Illustrate the rejection method
    [p,t]=distmesh2d_ill_reject_method(fd,@hscale_func,hmin,bbx,fpts,binf);    
end
%-------------------------------------------------------------------------------
figure
[p,t]=distmesh2d(fd,@hscale_func,hmin,bbx,fpts,binf);%not relates to p above. p is n x 3, t is triangulation of dimen m x 3

%-------------------------------------------------------------------------------
if dbg_flg == 1 
    hold on
    plot(fpts(:,1),fpts(:,2),'dg','linewidth',2,'markersize',3)
    lbl_fmt_fig('X','Y','')
    axis equal
    
end

%-------------------------------------------------------------------------------
% Construct the h-scaling function:
%   The scaling function is basic. For each input point we find the closest
%   Inner point and closest outer point. Then we construct a linear
%   function that changes from hin to hout based on where it is between the
%   two points. 
function hsc = hscale_func(p,binf)

%-------------------------------------------------------------------------------
% Extract the boundary data
hin   = binf.hin;
hout  = binf.hout;
houtv = binf.houtv;
pin   = binf.pin;
pout  = binf.pout;

%-------------------------------------------------------------------------------
% Find the closest points for each p
closest = zeros(size(p,1),2);
for n = 1:size(p,1)
    [~,closest(n,1)] = min( ( pin(:,1)-p(n,1)).^2 + ( pin(:,2)-p(n,2)).^2 );
    [~,closest(n,2)] = min( (pout(:,1)-p(n,1)).^2 + (pout(:,2)-p(n,2)).^2 );
end
hout = houtv(closest(n,2));

%-------------------------------------------------------------------------------
% Calculate the distance between the Inner to Outer closest points    
din2out = sqrt(sum( (pout(closest(:,2),:) - pin(closest(:,1),:)).^2,2));

%-------------------------------------------------------------------------------
% Construct a linear fit that scales between 1 and the ratio of the
% h-values
if hin < hout
    %---------------------------------------------------------------------------
    % Interior boundary spacing is smaller than Outer
    %---------------------------------------------------------------------------
    % Determine the distance of all ps to the closest interior points
    din2p   = sqrt(sum( ( p - pin(closest(:,1),:)).^2,2));
    
    %---------------------------------------------------------------------------
    % Calculate the scale for each: 
    % 
    %           scale = 1*(1-t) + hout/hin*t
    % 
    % where t is the distance from the inner boundary to p relative to the
    % outer boundary 
    t   = din2p ./ din2out;
    hsc = 1*(1-t) + hout/hin*t;    
        
else
    %---------------------------------------------------------------------------
    % Outer boundary spacing is smaller than Inner
    %---------------------------------------------------------------------------
    % Determine the distance of all ps to the closest outer points
    dout2p   = sqrt(sum( ( p - pin(closest(:,1),:)).^2,2));
    
    %---------------------------------------------------------------------------
    % Calculate the scale for each: 
    % 
    %           scale = 1*(1-t) + hout/hin*t
    % 
    % where t is the distance from the inner boundary to p relative to the
    % outer boundary 
    t   = dout2p ./ din2out;
    hsc = 1*(1-t) + hin/hout*t;
end

%-------------------------------------------------------------------------------
% Set a maximum stretch of the interior points
K       = binf.K;
IN      = find( inpolygon(p(:,1),p(:,2),pin(K,1),pin(K,2)) == 1);
% whos
% figure
% hold on
% plot(pin(:,1),pin(:,2),'.g','markersize',12)
% plot(pin(K,1),pin(K,2),'-r')
% plot(p(IN,1),p(IN,2),'.c','markersize',4)
if isempty(IN) == 0
    hsc(IN) = min( [hsc(IN) 20+0*IN],[],2);
end

% %-------------------------------------------------------------------------------
% figure
% plot3(p(:,1),p(:,2),hsc,'.k','markersize',8)
% clear