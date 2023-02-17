%-------------------------------------------------------------------------------
%
% Add the top and bottom of the triangulation
% 
%-------------------------------------------------------------------------------
function [t,p] = addtopbot_surftri(t,bp,bot_fps,top_fps,dbg_flg)

if nargin < 2
    error('stop')
elseif nargin < 3
    bot_fps = [];
    top_fps = [];
    dbg_flg = 0;
elseif nargin < 4
    top_fps = [];
    dbg_flg = 0;
elseif nargin < 5    
    dbg_flg = 0;
end

%-------------------------------------------------------------------------------
% Separate out the top and bottom edges
topi = find( abs( bp(:,3) - max(bp(:,3)) ) < 0.001);
boti = find( abs( bp(:,3) - min(bp(:,3)) ) < 0.001);

%-------------------------------------------------------------------------------
% Construct a triangulation of the top
if size(top_fps,1) == 0
    [ptop,ttop] = construct_distmesh_polygon_nofixed(bp(topi,1:2),0);
else
    [ptop,ttop] = construct_distmesh_polygon_innerfixed(bp(topi,1:2),top_fps,0);
end
ptop = [ptop max(bp(:,3))+0*ptop(:,1)];

%-------------------------------------------------------------------------------
% Construct a triangulation of the bottom
if size(bot_fps,1) == 0
    [pbot,tbot] = construct_distmesh_polygon_nofixed(bp(boti,1:2),0);
else
    [pbot,tbot] = construct_distmesh_polygon_innerfixed(bp(boti,1:2),bot_fps,0);
end
pbot = [pbot min(bp(:,3))+0*pbot(:,1)];

%-------------------------------------------------------------------------------
% Combine the triangles
% Combine the sides with the top
[p,t] = combine_triangulations(bp,ptop,t,ttop,0);
% Combine the sides and top with the bottom
if dbg_flg == 1
    [p,t] = combine_triangulations(p,pbot,t,tbot,1);
else
    [p,t] = combine_triangulations(p,pbot,t,tbot,1);
end

%-------------------------------------------------------------------------------    
% Plot the surface of the chest wrapped back up
if dbg_flg == 1    
    figure
    trisurf(t,p(:,1),p(:,2),p(:,3), 'Facecolor','cyan','FaceAlpha',0.8); axis equal;
    hold on    
    % saveas(gcf,'figs/step4_distmesh3Douter_added_top_bot','png')    
end