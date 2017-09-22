function  [Fx] = ForceX1_plate(gamma,Vf,Wf,rhoair,u,w,l,alpha,L, omega)

% Sphere static drag equations
% Calcul de la force selon X
% g=9.81;
% the drag is parallel to the relative velocity
% the lift is perpendical

% compute the relative velocity
vrel=sqrt((Vf-u).^2 + (Wf-w).^2);
% compute the relative Reynolds number
Cl = Cl_plate(gamma,gamma) + Clr_plate(vrel,omega);
Cd = Cd_plate(gamma,gamma) + Cdr_plate(vrel,omega);
% 
% k=(rhoair*Vf^2*l*L)/(2*g); %ici cest k=K/m, le m est dans les motion plate Shimoi
% Fx_Holmes1= k*g*((1-(u/Vf)).^2 + (w/Vf).^2)*(Cd.*cosd(-alpha)-Cl*sind(-alpha)); 

Fx= (0.5*rhoair*l*L*(vrel.^2.0)*(Cd.*cosd(alpha)-Cl*sind(alpha)));

end