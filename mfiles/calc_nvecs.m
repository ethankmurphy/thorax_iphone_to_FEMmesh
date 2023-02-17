function [ns,a,b] = calc_nvecs(varargin)

if nargin == 1
    p = varargin{1}.node;
    t = varargin{1}.face;
elseif nargin == 2
    p = varargin{1};
    t = varargin{2};
end


%---------------------------------------------------------------------------
a = p(t(:,2),:) - p(t(:,1),:);
b = p(t(:,3),:) - p(t(:,1),:);
% the normal is the cross product of two vectors
% a x b = (a2*b3-a3*b2)i + (a3*b1-a1*b3)j + (a1*b2-a2*b1)k
ns = [ ...
    (a(:,2).*b(:,3) - a(:,3).*b(:,2)) ...
    (a(:,3).*b(:,1) - a(:,1).*b(:,3)) ...
    (a(:,1).*b(:,2) - a(:,2).*b(:,1))];
% Normalize
nlens = sqrt(sum(ns.^2,2));
ns = ns ./ repmat(nlens,1,3);
% Normalize the a and b
alens = sqrt(sum(a.^2,2));
blens = sqrt(sum(b.^2,2));
a     = a ./ repmat(alens,1,3);
b     = b ./ repmat(blens,1,3);
