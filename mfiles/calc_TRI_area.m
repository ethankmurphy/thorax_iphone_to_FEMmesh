function areaTRI = calc_TRI_area(TRI,nodes,useparfor)

if nargin == 2
    useparfor = 0;
end

%-------------------------------------------------------------------------------
% Find the overall area based on the dimensions and assumptions its a
% circle
% rad      = abs(min(nodes(:)));
% tot_area = pi*rad^2

if (size(nodes,2) == 2) || (abs(std(nodes(:,3))) < 1e-6)
    %-------------------------------------------------------------------------------
    % Loop through each element and calculate the area
    % The area should be 1/2 * abs(det( two vectors ))
    areaTRI = zeros(size(TRI,1),1);
    if useparfor == 1
        parfor n = 1:size(TRI,1)
            v1   = [nodes(TRI(n,1),1:2)'; 0];
            v2   = [nodes(TRI(n,2),1:2)'; 0];
            v3   = [nodes(TRI(n,3),1:2)'; 0];
            areaTRI(n) = norm(cross( (v2-v1), (v3-v1) ));
        end
    else
        for n = 1:size(TRI,1)
            v1   = [nodes(TRI(n,1),1:2)'; 0];
            v2   = [nodes(TRI(n,2),1:2)'; 0];
            v3   = [nodes(TRI(n,3),1:2)'; 0];
            areaTRI(n) = norm(cross( (v2-v1), (v3-v1) ));
        end
    end
    areaTRI = areaTRI/2;
    
elseif size(nodes,2) == 3
    %-------------------------------------------------------------------------------
    % Loop through each element and calculate the area
    % The area should be 1/2 * abs(det( two vectors ))
    areaTRI = zeros(size(TRI,1),1);
    if useparfor == 1
        parfor n = 1:size(TRI,1)
            v1   = nodes(TRI(n,1),:)';
            v2   = nodes(TRI(n,2),:)';
            v3   = nodes(TRI(n,3),:)';
            areaTRI(n) = norm(cross( (v2-v1), (v3-v1) ));
        end
    else
        for n = 1:size(TRI,1)
            v1   = nodes(TRI(n,1),:)';
            v2   = nodes(TRI(n,2),:)';
            v3   = nodes(TRI(n,3),:)';
            areaTRI(n) = norm(cross( (v2-v1), (v3-v1) ));
        end
    end
    areaTRI = areaTRI/2;
end