%-------------------------------------------------------------------------------
%
% Get the centers of a mesh. 
%  
% Usage: get_tcs(p,t) or get_tcs(msh)
% 
%     * where p and t are node and triangle matrices assumed size Nnx3 and
%       Ntx3 where Nn are the number of nodes and Nt is the number of
%       triangles.
%     * where msh is a structure array with two fields .node and .face that
%       are nodes and elements/triangles defined the same as p and t.
%
%-------------------------------------------------------------------------------
function tcs = get_tcs(varargin)

if nargin == 1
    p = varargin{1}.node;
    t = varargin{1}.face;
elseif nargin == 2
    p = varargin{1};
    t = varargin{2};
end
    

tcs = 1/3*(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:));
