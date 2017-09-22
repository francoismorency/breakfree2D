function  [Cl] = Cl_plate(gamma,theta)
%%
% lift force equations from normal force
% Calculation of the lift coefficient in term of the angle of attack and
% the normal coefficient
global dim

if (dim == 2)
% version 2D infinite span plate
CN_theta = CN_2D(gamma);
else
%version 3D square plate
CN_theta = CN_plate(gamma);
end

Cl=CN_theta.*cosd(mod(theta,180)); %%(eq.40b)

% a='Cl_Homes'

end