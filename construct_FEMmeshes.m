%--------------------------------------------------------------------------
%
%     construct_run_gmsh
%
%--------------------------------------------------------------------------
clear
clc
close all

%--------------------------------------------------------------------------
% Set parameters
hbnd      = 5;    % h-value in background (mm)
hel       = 0.5;    % h-value on electrodes (mm)
delec     = 10;     % Diameter (mm)
repo_pth  = [];      % Path of where processed data will be saved/loaded
dbg_flg   = 1;

%--------------------------------------------------------------------------
% Set the path
addpath(genpath('mfiles'))
%---------------------------------------------
% Required software
%   1. gmsh - below is a string to the executable, notice the space after
%             the gmsh in the string
% rungmsh_str = '/jumbo/digihisto/Ethan/software/gmsh-4.3.0-Linux64/bin/gmsh ';
rungmsh_str = 'C:\gmsh-4.8.4-Windows64\gmsh ';
%---------------------------------------------
%   2. distmesh - This is a matlab simple (mainly) 2D meshing software. 
%                 You just need to download it, https://github.com/ionhandshaker/distmesh
%                 and then add the path to it. 
addpath S:\digihisto\Ethan\software\distmesh

%--------------------------------------------------------------------------
% Set the scan name and path to be analyzed. This is of a very particular
% format. It is a processed 3D iPhone scan from the Heges app using the
% Github repository: https://github.com/ethankmurphy/EIT_3D_thorax_scans
if dbg_flg == 0
    [file,path] = uigetfile('*.MAT','Pick a Mat-file of a processed iPhone 3D Thorax Scan (suffix of ''_boundmesh_elecs'')');   % Select the file of interest
else    % Debugging
    path = 'example';
    file = '1675357621o1471639_boundmesh_elecs';
end
if isempty(repo_pth)
    repo_pth = path;
end
addpath mfiles


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Load the boundary information and the electrodes
eval(['load ',path,'/',file,' bmsh elecs colobj tri bndfit' , ...
    ' total_time order_fit_time click_time load_crop_time'])
tri0    = tri;
%-------------------------------------
% Convert to units of mm
bmsh.node = bmsh.node*1000;
elecs     = elecs*1000;
locs      = colobj.Location*1000;
colobj    = pointCloud(locs, 'Color', colobj.Color );
%------------------------------------
if dbg_flg == 1
    figure; hold on
    plot_colobj_tri(colobj,tri)
    trisurf(bmsh.tri,bmsh.node(:,1),bmsh.node(:,2),bmsh.node(:,3),'facecolor','cyan','linestyle','none', ...
        'facealpha',0.5)
    plot3(1.01*elecs(:,1),1.01*elecs(:,2),1.01*elecs(:,3),'.m','markersize',20)
    axis equal
    for n = 1:size(elecs,1)
        text(1.03*elecs(n,1),1.03*elecs(n,2),1.03*elecs(n,3), ...
            num2str(n), ...
            'color','yellow')
    end
    title('example scan')
    view([10 6])
end

%----------------------------------------------------------------------
% Construct or load the outer surface triangulation with encoded
% electrodes
tic
ply_matnam = [file(1:end-16),'_encoded_outersurf'];
if filechecker(path,[ply_matnam,'.mat']) == 0
    %------------------------------------------------------------------
    % Construct 3D electrode boundary points tangent to the surface of
    % the boundary mesh
    [el3D_bndpts,elns] = add3D_elec_bndpts_v2(elecs,delec,bmsh,0);

    %------------------------------------------------------------------
    [t,bp,el3D_bndpts,area_errs] = unwrap_bnd_v3(el3D_bndpts,delec,bmsh,elecs,bndfit,elns,0);

    %------------------------------------------------------------------
    % Add a top and bottom surface and check that the surface is water
    % tight
    [t,bp] = addtopbot_surftri(t,bp);

    %------------------------------------------------------------------
    % Save the data
    encode_outersurf_time = toc
    eval(['save ',repo_pth,'/',ply_matnam,' t bp el3D_bndpts area_errs encode_outersurf_time'])
else
    %----------------------------------------------------------------------
    % Load the data
    eval(['load ',repo_pth,'/',ply_matnam,' t bp el3D_bndpts area_errs encode_outersurf_time'])
end


%-------------------------------------
if dbg_flg == 1
    figure; hold on
    plot_colobj_tri(colobj,tri,0.5)
    trisurf(t,bp(:,1),bp(:,2),bp(:,3),'facecolor','cyan', ... 'linestyle','none', ...
        'facealpha',0.5)
    plot3(1.01*elecs(:,1),1.01*elecs(:,2),1.01*elecs(:,3),'.m','markersize',20)
    axis equal
    for n = 1:size(elecs,1)
        text(1.03*elecs(n,1),1.03*elecs(n,2),1.03*elecs(n,3), ...
            num2str(n), ...
            'color','yellow')
    end   
    view([10 6])
    drawnow;pause(1)
end


%----------------------------------------------------------------------
%----------------------------------------------------------------------
% Construct the 3D FEM mesh
tic
ply_matnam = [file(1:end-16),'_3D_thrx'];
if filechecker(path,[ply_matnam,'.mat']) == 0
    %------------------------------------------------------------------
    [msh,eareas] = construct_FEMmsh_func(bp,t,el3D_bndpts,delec,hbnd,hel,elecs,rungmsh_str,0);

    %------------------------------------------------------------------
    % Save the data
    FEMmsh_time = toc
    eval(['save ',repo_pth,'/',ply_matnam,' msh eareas FEMmsh_time'])
else
    %------------------------------------------------------------------
    % Load the data
    eval(['load ',repo_pth,'/',ply_matnam,' msh eareas FEMmsh_time'])
end

%----------------------------------------------------------------------
% Plot the mesh with its interior surfaces
figure;set_fig_relsiz(0.5);hold on
simpplot(msh.node*1000,msh.elem,'p(:,2)>0|p(:,1)>0');hold on
%--------------------------------------------------------
locs        = colobj.Location;
ris         = find( (colobj.Location(:,1) < 0 ) & (colobj.Location(:,2) < 0 ));
[locs,tri]  = update_tris_due_to_rmvnodes(locs,tri0,ris);
locs(:,1:2) = locs(:,1:2)*1.02;
cols        = colobj.Color;
cols(ris,:) = [];
colobjtmp   = pointCloud(locs, 'Color', cols);
plot_colobj_tri(colobjtmp,tri,1)
%--------------------------------------------------------
for n = 1:length(msh.elec)
    ftmp   = msh.face(msh.elec{n},:);
    trisurf(ftmp,1.05*msh.node(:,1)*1000,1.05*msh.node(:,2)*1000,msh.node(:,3)*1000, ...
        'FaceColor','red','FaceAlpha', 1,'linestyle','none');
end
view(3)
%--------------------------------------------------------
axis equal;axis off
view([-16 8]);camlight right
saveas(gcf,'example_processed_FEMmesh','png')

%----------------------------------------------------------------------
disp('*******************************************************')
disp(['    Load/Crop time: ',num2str(round(load_crop_time,1)),' s'])
disp(['   Elec Click time: ',num2str(round(click_time,1)),' s'])
disp(['   Elec Order time: ',num2str(round(order_fit_time,1)),' s'])
disp(['Outer Trisurf time: ',num2str(round(encode_outersurf_time,1)),' s'])
disp(['     FEM mesh time: ',num2str(round(FEMmsh_time,1)),' s'])
disp(['        Total time: ',num2str(round( ...
    load_crop_time + click_time + order_fit_time + encode_outersurf_time + FEMmsh_time,1)),' s'])
disp(' ')

