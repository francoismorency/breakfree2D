function  [CM] = Moment_Coefficient(gamma)


% Calcul du coefficient de moment par rapport � l'angle d'attaque
% gamma
% LxHolmes= position adimensionn�e du centre de pression par rapport au
% centre de la plaque.
% CN_theta= Normal coefficient (gamma)
% gamma= angle of attack

Lx = Centre_of_Pressure_plate(gamma);
CN_theta = CN_plate(gamma);

CM=CN_theta.*Lx; %%(eq.41)


end