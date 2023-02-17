%-------------------------------------------------------------------------------
%
%                            read_gmsh_mesh
%
%-------------------------------------------------------------------------------
%
% 2.  Extract the nodes and elements from the gmsh mesh file
%
%-------------------------------------------------------------------------------
function msh =                ...
    read_gmsh_mesh_v3(           ...
    mshfname                , ...
    dbg_flg)


%-------------------------------------------------------------------------------
% Open the mesh file
fid = fopen([mshfname,'.msh']);

tline    = fgets(fid);
stoploop = 0;
while (ischar(tline)) && (stoploop == 0)
    
    %---------------------------------------------------------------------------
    % Get for the nodes
    if isempty( strfind(tline,'$Nodes') ) == 0
        %-----------------------------------------------------------------------
        % The next line is the number of nodes
        numnode = str2num( fgets(fid) );
        
        %-----------------------------------------------------------------------
        % Read in every node
        tmp  = textscan(fid,'%f',4*numnode);
        node = reshape(tmp{1},4,numnode)';
        node = node(:,2:end);
    end
    
    %---------------------------------------------------------------------------
    % Get for the Elements
    if isempty( strfind(tline,'$Elements') ) == 0
        %-----------------------------------------------------------------------
        % The next line is the number of elements
        numelem = str2num( fgets(fid) );
        
        %-----------------------------------------------------------------------
        % Read in every element, but stop when we reach the first
        % tetrahedral element.
        face  = zeros(numelem,3);
        q     = 1;
        eltyp = 0;
        while eltyp ~= 4
            %-------------------------------------------------------------------
            % Get the line and convert it to a numeric vector
            tmpnum  = str2num( fgets(fid) );
            tmpelem = tmpnum(1:end);  
            try
            eltyp   = tmpelem(2);
            catch
                rethrow(eltyp);
            end
            numelem = numelem-1;            
            %-------------------------------------------------------------------
            % Record the triangles
            if eltyp == 2
                face(q,:)  = tmpelem(6:8);
                q          = q+1;
            end
        end
        
        %-----------------------------------------------------------------------
        % Read the rest of the tetrahedral elements with textscan        
        tmp    = textscan(fid,'%f',9*numelem);                
        tmpels = reshape(tmp{1},9,numelem)';        
        elem   = [ ...
            tmpelem(6:9); ...
            tmpels(:,6:9)];
          
        stoploop = 1;
    end
    
    tline = fgets(fid);
end
fclose(fid);

% %------------------------------------------------------------------------
% % Only take the tetrahedron elements
% inds = find( elem(:,1) == 4);
% elem = elem(inds,5:8);

%--------------------------------------------------------------------------
% Construct the msh structure array
msh.node = node;
msh.elem = elem;
msh.face = face(1:(q-1),:);

%--------------------------------------------------------------------------
% Debugging
if dbg_flg == 1
    figure
    plot3(node(:,1),node(:,2),node(:,3),'ok')    
end
end