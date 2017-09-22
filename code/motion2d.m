function [dx] = motion2d (t,x);

%%%%% Fonction basee sur la fonction originale 'motion_plate_Shimoi'
%%%%% t time, aerodynamic coefficients do not depend on time, but required by ode
%%%%% x vector of unknown

% fonctions requises pour le fonctionnement du code : 

%alpha, CN_Holmes, Clift, Clr_holmes, Cd_Holmes, Cdr_Holmes, Centre_of_Pressure_Holmes, Moment_Coefficient,
%Moment Holmes, Cmr_Holmes, ForceX1_Holmes, ForceZ1_Holmes, figXZ1.


% La resolution de ce programme repose sur des corrï¿½lations donnant un
% coefficient de force normal Cn et la position du centre de position du
% centre de pression c, etablies par J.D Holmes dans "investigation of
% plate-type windborne debris --PartII : Computed trajectories, 2006".

%Hypotheses : 
%L'ecoulement est uniforme et de vitesse constante Vf, la plauque ne
%perturbe pas cet ecoulement.
%La resolution a lieu dans un champ de gravite constant, g=9.81.

% derivative and variables
% x(1) = x
% x(2) = Vx=u
% x(3) = z
% x(4) = Vz=w
% x(5) = theta, the angle of the geometry in x,z coordinate
% x(6) = omega = dtheta/dt

%%this is global variable
global l L e dim Vf Wf muair rhoair g rhoice dsp;

% definition of the freestream variables
%Vf = 70.0 ; % freestream  x velocity
%Wf = 0.0 ; % freestream w velocity
%muair=1.e-5; % viscosité de l'air
%rhoair=1.24;  %% densitï¿½ de l'air (milieu dans lequel les calculs sont faits) (kg/m^3)
%g=9.81;      %% gravity
%rhoice=57.25*0.45359237/(0.3048^3); %% 917.0570 masse vol de la glace kg/m^3 annoncï¿½ par Koji Shimoi
%rhoPlastique=1120; %% masse volumique du plastique(kg/m^3)
% defined in TraceXZ, should be read in a file eventually
%Vf=1.0;
%Wf=0.0;
%muair=1.14e-4;
%rhoair=1.0;
%g= 2.25e-4;%3.75e-4 ou 2.25e-4;




%%Constantes pour une plaque plane
% defined in traceXZ
%dim=2;  %dim=2 2D Cn values, else 3D Cn values
%l=6*0.0254 ;     %% corde de la plaque (m)`
%L=6*0.0254 ;     %% envergure de la plaque (m)`
%e=0.04*0.0254 ;    %% ï¿½paisseur de la plaque (m)
%rhoice=660;%2.38e2
%l=1.0 ;     %% corde de la plaque (m)`
%L=10.0 ;     %% envergure de la plaque (m)`
%e=0.14966 ;    %% ï¿½paisseur de la plaque (m)


%%constante pour une sphère
%dsp=0.04 % diameter
%uts=g*dsp^2.*(rhoice-rhoair)/(18.*muair) %terminal velocity
% the terminal velocity is not needed, but is usefull 
% the terminal velocity assume stoke flow
% see equation 3-18, "Bubbles, Drops, and Particles" 
% Clift, Grace, Weber

%Conditions de tests :

%Tachikawa : plaque carree de plastique rho=1120 kg/m^3, l=L=0.04m,
%e=0.002m,Vf=9.18m.s-1, K=rhoair*(l*Vf)^2/(2*m*g)=2.3 => rhoair = 1.195,
%theta 0 = 15 et 45 degres, a comparer sur X=1m

%Shimoi SFP : rhoplaque=rhos Shimoi, l=L=6*0.0254, e=0.04*0.0254,
%Vf=71.5264, theta=0, theta=90

% Shimoi RFP12 : rhoplaque=rhos Shimoi, l=0.5*L=6*0.0254, e=0.04*0.0254,
%Vf=71.5264, theta=0, theta=-10

% Shimoi RFP6 : rhoplaque=rhosice, 0.5*l=L=6*0.0254, e=0.04*0.0254,
%Vf=71.5264, theta=0



% Masse de la plaque (kg)
mice=rhoice*l*L*e;
% Masse de la sphere
%mice=rhoice*4./3.*pi*(dsp/2.)^3.;

% moment d'inertie pour une plaque rectangulaire de corde l , envergure L et d'ï¿½paisseur e. ((kg.m^2) par rapport ï¿½ l'axe Y
% $$ I = m_{ice} \frac{(e^2+l^2)}{12}
I=mice*(e^2+l^2)/12.0;
% moment d'inertie pour une sphere
%I=mice*2./5.*(dsp/2.)^2.;


% angle du vent relatif relatif
alpha1=alpha(x(2), x(4), Vf, Wf);

% angle d'attaque total = alpha+theta
gamma =  mod(x(5)+ alpha1,180) ;


%le x(6) sert ï¿½ prendre en compte l'effet Magnus, par dï¿½fault l'effet
%Magnus n'est pas implï¿½menter (mis en commentaire). Pour le prendre en
%compte ajouter les coefficient Clr_Holmes, Cdr_Holmes, Cmr_Holmes dans les
%fonctions ForceX1_Holmes, ForceZ1_Holmes, Moment_Holmes. (ces coefficient
%sont dï¿½jï¿½ prï¿½sents en commentaires dans ces fonctions).

% force en X plaque plane
if (dim > 1)
    FD =ForceX1_plate(gamma,Vf,Wf,rhoair,x(2),x(4),l,alpha1,L,x(6));
% force en X, sphere
else
    FD =ForceX1_sphere(gamma,Vf,Wf,rhoair,muair,x(2),x(4),dsp,alpha1,x(6));
end

% force en Z, plaque plane
if (dim > 1)
    FP = ForceZ1_plate(gamma,Vf,Wf,rhoair,x(2),x(4),l,alpha1,L,x(6));
% force en Z, sphere
else
    FP = ForceZ1_sphere(gamma,Vf,Wf,rhoair,muair,x(2),x(4),dsp,alpha1,x(6));
end

% Moment, plaque plane
if (dim > 1)
    Mt = Moment_plate(gamma ,Vf,Wf, rhoair,x(2),x(4),l,L,x(6),I,mice);
% Moment, sphere
else
    Mt = 0.0d0;
end

% Trajectoire X
dx(1,1)=x(2);       %dx/dt=u     
dx(2,1)=FD/mice;      %du/dt=Fx/mice

% trajectoire Z
dx(3,1)=x(4);      %dz/dt=w       
dx(4,1)=FP/mice - g;   %dw/dt=Fz/mice - g


% theta angle avec l'horizontal
dx(5,1)=x(6)*180/pi; %dtheta/dt=omega (theta en degrï¿½)     
dx(6,1)=(Mt/(I)); %% domega/dt=My/Iy




end



