function  [Cd] = Cd_plate(gamma,theta)

% Holmes equations
% Calculation of the drag coefficient in term of the angle of attack and
% the normal coefficient
global dim

if (dim == 2)
%version 2D
CN_theta = CN_2D(gamma);
else
CN_theta = CN_plate(gamma);
end


Cd=0.1+ CN_theta'.*sind(mod(theta,180)); %%(eq.40a)


end