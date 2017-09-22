function  [Mt] = Moment_plate(gamma,Vf,Wf,rhoair,u,w,l,L,omega,I,mice)

%%Holmes  equations
%%Calcul du moment autour de l'axe Y (moment de tangage)
%%gamma= angle d'attaque total
% g=9.81;
% K=(rhoair*Vf^2*l^2)/(2*mice*g);
% Ln=g*l/Vf^2;
% In=I/(mice*l^2);
% U=g*u/Vf^2;
% W=g*w/Vf^2;
% compute the relative velocity
vrel=sqrt((Vf-u).^2 + (Wf-w).^2);


%Pour prendre en compte l'effet Magnus rajouter Cmr
CM =Moment_Coefficient(gamma) + Cmr_plate(vrel,omega,u,w);

% Mt_Holmes=(K/(Ln*In))*((1-U)^2+W^2)*CM_Holmes;
Mt= 0.5*(CM*rhoair*((l*L)*l)*vrel.^2);



end