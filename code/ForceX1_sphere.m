function  [Fx_sp] = ForceX1_sphere(gamma,Vf,Wf,rhoair,muair,u,w,dsp,alpha, omega)

% Sphere static drag equations
% Calcul de la force selon X
% g=9.81;
% the drag is parallel to the relative velocity
% the lift is perpendical

% compute relative velocity, 
% function out = CD_sp(Re)
% Calcul du coefficient de trainée de la sphère
% the coherence of sign bewteen w and g must still be verified
vrel=sqrt((Vf-u).^2 + (Wf-w).^2);
% compute the relative Reynolds number
Re=rhoair*dsp*vrel/muair;

% compute de drag coefficient and lift
Cd = CD_sp(Re);
Cl = 0.0d0;


% compute the force in absolute x direction
% for generality, wind relatie angle is taken into account 
% alpha is the wind relative angle
Afront= pi*dsp.^2/4.0;
Fx_sp= 0.5*rhoair*Afront*vrel.^2.*(Cd.*cosd(alpha)-Cl*sind(alpha));

end