function [p,t] = combine_triangulations(p1,p2,t1,t2,check_freeb,roundflg)

if nargin < 5
    % If the final variable, check_freeb, was not input to the function,
    % then assign it a value of zero to skip the free boundary check.
    check_freeb = 0;
    roundflg    = 0;
elseif nargin < 6
    roundflg    = 0;
end

if roundflg > 0
    p1 = round(p1,roundflg);
    p2 = round(p2,roundflg);
end



%p-points %t-triangulation
%-------------------------------------------------------------------------------
% Combine the triangles
np1 = size(p1,1);
np2 = size(p2,1);
t   = [t1; (t2+np1)];
p   = [p1; p2];

%-------------------------------------------------------------------------------
% Connect the mesh 2 (p2,t2) with mesh 1 (p1,t1). It is assumed that they in
% fact intersect nicely, but just have repeated edges.
for n = 1:np2
    dists     = sqrt( sum( (p1 - repmat(p2(n,:),np1,1)).^2,2));
    [mdist,i] = min(dists);
    if mdist < 1e-3
        % Replace the triangles with intersecting nodes to correspond to the
        % nodes from mesh 1
        for k = 1:size(t,2)
            i0      = find( t(:,k) == (n+np1));
            t(i0,k) = i;
        end
    end
end

%-------------------------------------------------------------------------------
% Remove unused nodes
if size(t,2) == 3
    [used_nds,tmp,uis] = unique([t(:,1); t(:,2); t(:,3)]);
elseif size(t,2) == 4
    [used_nds,tmp,uis] = unique([t(:,1); t(:,2); t(:,3); t(:,4)]);
end
p = p(used_nds,:);

%-------------------------------------------------------------------------------
% Relabel the mesh elements (second method)
ne         = size(t,1);
elnew      = t;
if size(t,2) == 3
    elnew(:,1) = uis(       1:ne);
    elnew(:,2) = uis((  ne+1):2*ne);
    elnew(:,3) = uis((2*ne+1):3*ne);
elseif size(t,2) == 4
    elnew(:,1) = uis(       1:ne);
    elnew(:,2) = uis((  ne+1):2*ne);
    elnew(:,3) = uis((2*ne+1):3*ne);
    elnew(:,4) = uis((3*ne+1):4*ne);
end
t          = elnew;

%----------------------------------------------------------------------
% Look for and remove degenerate triangles
tareas = calc_TRI_area(t,p);
idegen = find(tareas < 1e-10);
t(idegen,:) = [];


%--------------------------------------------------------------------------
% Check that the freeboundary is empty (if requested)
if check_freeb == 1
    DT = triangulation(t, p);
    ks = freeBoundary(DT);
    %----------------------------------------------------------------------
    if size(ks,1) > 0
        %------------------------------------------------------------------
        % Try to fix it by first removing degenerate triangles
        areaTRI = calc_TRI_area(t,p,0);
        figure
        plot(sort(areaTRI))
        set(gca,'yscale','log')
        %-----------------------------
        isrmv      = find(areaTRI < 1e-8);
        t(isrmv,:) = [];
        %------------------------------------------------------------------
        % Remake the freeBoundary
        DT = triangulation(t, p);
        ks = freeBoundary(DT);

        if size(ks,1) > 0
            %-----------------------------
            % Ok, can't be easily fixed :(, report it
            size(ks,1)
            sort(ks,2)
            p(unique(ks(:)),:)
            
            abs(sqrt(comp_pairwise_distmat(p(unique(ks(:)),:))))
            
            figure
            trisurf(ks,p(:,1),p(:,2),p(:,3))

            figure
            hold on
            trisurf(t1,p1(:,1),p1(:,2),p1(:,3),'facecolor','cyan','facealpha',0.2)
            trisurf(t2,p2(:,1),p2(:,2),p2(:,3),'facecolor','green','facealpha',0.2)
            trisurf(ks,p(:,1),p(:,2),p(:,3),'facecolor','red')
            error('non empty free boundary, a.k.a. problem!')
        else
            disp('Surface Triangulation is water tight')
        end
    else
        disp('Surface Triangulation is water tight')
    end

end
