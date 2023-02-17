function [uv1,uv2] = get_two_orthvecs(nvec,alignvec)


%-------------------------------------------------------------------------------
% Unit vector 1
vv1 = rand(1,3);
uv1 = cross(nvec,vv1);
uv1 = uv1 / norm(uv1);
%-------------------------------------------------------------------------------
% Unit vector 2
uv2 = cross(nvec,uv1);
uv2 = uv2 / norm(uv2);

%-------------------------------------------------------------------------------
% If there is also a specified alignment vector, then do a small
% optimization to best-align the uv1 vector with this vector
if nargin == 2
    %-------------------------------------------------------------------------------
    % Rotate the uv1 vector so its most aligned with the x-axis
    bestt         = fminbnd( @(t) rot_uv1basis_align(t,uv1,uv2,alignvec),0,2*pi);
    [tmp,uv1,uv2] = rot_uv1basis_align(bestt,uv1,uv2,alignvec);
    
end



%-------------------------------------------------------------------------------
%
% rot_uv1basis_align
%
%-------------------------------------------------------------------------------
function [delt,uv1r,uv2r] = rot_uv1basis_align(t,uv1,uv2,alignvec)

%-------------------------------------------------------------------------------
% Rotate a vector representing 100% in the uv1 direction to be some angle
% off from this
R   = [cos(t) -sin(t); sin(t) cos(t)];
vr1 = R*[1;0];
vr2 = R*[0;1];

%-------------------------------------------------------------------------------
% Get the rotated vectors
uv1r = vr1(1) * uv1 + vr1(2)*uv2;
uv2r = vr2(1) * uv1 + vr2(2)*uv2;
% [norm(uv1r) norm(uv2r) dot(uv1r,uv2r)]
%-------------------------------------------------------------------------------
% Calculate the angle between uv1r and [1 0 0];
% delt = abs( acos(dot(uv1r,[1 0 0])) );
delt = abs( acos(dot(uv1r,alignvec)) );