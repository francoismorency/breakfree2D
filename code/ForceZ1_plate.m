function  [Fz] = ForceZ1_plate(gamma,Vf,Wf,rhoair,u,w,l,alpha,L,omega) ;
%%
% Holmes  equations for forces
% Calcul de la force selon Z
% w= dz/dt

% compute the relative velocity
vrel=sqrt((Vf-u).^2 + (Wf-w).^2);

%Pour prendre en compte l'effet Magnus rajouter Clr et Cdr
 Cl = Cl_plate(gamma,gamma) + Clr_plate(vrel,omega);
 Cd = Cd_plate(gamma,gamma) + Cdr_plate(vrel,omega);
 

 Fz=(0.5*rhoair*l*L*(vrel.^2)*(Cl.*cosd(alpha)+Cd.*sind(alpha)));



end

