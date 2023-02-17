%--------------------------------------------------------------------------
%
% Evaluate the distance between points (p) and test points (testps).
% Specifically, for each test point we want to find the closest point (p)
% and its corresponding index.
% 
% One can obviously do this in a brute force way (methflg=1) where you loop
% through each test point and calculate a Euclidean distance, but when the
% sets of points get large this gets very slow. 
%
%
%
%--------------------------------------------------------------------------
function [ds,is] = eval_closest_pdist(p,testps,parflg,methflg,dbg_flg)



%--------------------------------------------------------------------------
if methflg == 1
    %----------------------------------------------------------------------
    % Brute force

    %----------------------------------------------------------------------
    ds = zeros(size(testps,1),1);
    is = zeros(size(testps,1),1);
    if parflg == 1
        parfor n = 1:size(testps,1)
            %-------------------------------------------------------------------
            % Find the closest tangent plane center to the current point
            [ds(n),is(n)] = min(sqrt( (p(:,1)-testps(n,1)).^2 + (p(:,2)-testps(n,2)).^2 + (p(:,3)-testps(n,3)).^2));

        end
    else
        if size(p,2) == 2
            for n = 1:size(testps,1)
                %-------------------------------------------------------------------
                % Find the closest tangent plane center to the current point
                [ds(n),is(n)] = min(sqrt( (p(:,1)-testps(n,1)).^2 + (p(:,2)-testps(n,2)).^2));
            end
        elseif size(p,2) == 3
            for n = 1:size(testps,1)
                %-------------------------------------------------------------------
                % Find the closest tangent plane center to the current point
                [ds(n),is(n)] = min(sqrt( (p(:,1)-testps(n,1)).^2 + (p(:,2)-testps(n,2)).^2 + (p(:,3)-testps(n,3)).^2));
            end
        end
    end
end