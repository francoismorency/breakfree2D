function [t,x] = TraceXZ
%% Plot trajectories in $XZ$ plan 
% Trace dans le plan $XZ$ la trajectorie d'un objet dans un ecoulement
% uniforme de vitesse $V_f$, $W_f$.
clear;

%% Initial conditions
%%
% * X0 position initiale en X
% * u0 vitesse initiale en X (positive de gauche a droite)
% * Z0 Poisiton initiale en Z
% * w0 Vitesse initiale en Z (positive vers le haut)
% * theta0 angle initial par rapport a l'horizontal de la plaque (cet angle
% est pris modulo 180)
% * omega0(vitesse angulaire initiale)


%%
%
% $$X_{in} = (X_0, u_0, Z_0, w_0, \theta_0, \omega_0 )$$
%
Xin=[0, 0, 0, 0, 60, 0];

%% simulation time and sampling time definition
%%
% * rotation positive => effet magnus diminuant la portance
% * rotation negative => effet magnus augmentant la portance
% * angle positif => prise au vent porteuse rotation positive
% * angle negatif => rotation negative
%%
% temps de simulation (s)
tmax=30.0;

%%
% note le pas de temps est variable dans ode23
% mais nous forçons l'écriture à un interval deltat (s)
% pas de temps
deltat=0.01;
%%
% vecteur du temps de simulation: sampling time
xspan=0.0:deltat:tmax; 

%% define flow and plate properties
% This is global variable needed for calculation that
% should be put in an input file eventually:
%%
% 
% * $V_f$ horizontal velocity;
% * $W_f$ vertical velocity (g direction);
% * $\mu_{air}$ air viscosity for sphere/cylinder;
% * $g$  gravity magnitude;
% * $\rho_{Shimoi}$ plate density;
% * dim either 2D or 3D (square plate) aerodynmic coefficients
% 
global Vf Wf l L muair rhoair g rhoice dim e dsp;

Vf=1.0;
Wf=0.0;
muair=1.14e-4;
rhoair=1.0;
g= 2.25e-4;%3.75e-4 ou 2.25e-4;
rhoice=660;%2.38e2
dim=3;  %dim=2 2D Cn values, dim=3 3D Cn values, dim=1 Cylinder value

%%
% * l : corde de la plaque (m)`
% * L : envergure de la plaque (m)`
% * e : epaisseur de la plaque (m)
l=1.0 ;     %% corde de la plaque (m)`
L=10.0 ;     %% envergure de la plaque (m)`
e=0.14966 ;    %% epaisseur de la plaque (m)

%%
% * dsp : diametre de la sphere
dsp=0.04; % diameter (m)

switch dim
    case 1
    disp('sphere');
    case 2
    disp('plate 2D correlation');
    case 3
    disp('plate 3D correlation');
    otherwise
    error('dim value should be 1,2 or 3. Actual dim=%i\n',dim)       
end


%% ODE numerical solution:
% compute object trajectory by solving ODE

options = odeset('OutputFcn',@odeplot);
[t, x]=ode23('motion2d', xspan, Xin,options);

%% Graphical outupt
% Mise des resultats sous forme graphique

% trajectoire dans le plan XZ

subplot(3,2,1), plot(x(:,1),x(:,3));
title('position');
%   axis([0 1 -0.21 0.35]);
xlabel('X (m)')
ylabel('Z (m)')
grid on;

% Vitesse selon X
subplot(3,2,2), plot(t,x(:,2));
title('u velocity m/s');
xlabel('t (s)')
ylabel('U (m/s)')
grid on;

% Vitesse angulaire
subplot(3,2,3), plot( t,(x(:,6)));
title('omega ');
xlabel('t(s)')
ylabel('omega (rad/s)')
grid on;

% angle d'attaque
 
%alpha1=alpha(x(:,2), x(:,4), Vf, Wf);
%gamma =  mod(x(:,5)+ alpha1,180); 
%FP = ForceZ1_Holmes(gamma(1),Vf,Wf,1.24,x(1,2),x(1,4),l,alpha1(1),L,x(1,6))
%FD =ForceX1_Holmes(gamma(1),Vf,Wf,1.24,x(1,2),x(1,4),l,alpha1(1),L,x(1,6))
%subplot(3,2,3), plot( t,gamma);
%title('AOA ');
%xlabel('t(s)')
%ylabel('deg')
%grid on;

% Evolution de l'angle de la plaque avec l'horizontal theta
subplot(3,2,4), plot( t,x(:,5));
title('orientation angle theta ï¿½')
xlabel('t (s)')
ylabel('theta (ï¿½)')
grid on;

% Evolution de l'angle alpha 
subplot(3,2,5), plot(t,alpha(x(:,2),x(:,4),Vf,Wf))
title('angle relatif alpha');
xlabel('t (x)')
ylabel('alphaï¿½')
grid on;

% vitesse en Z
subplot(3,2,6),
plot(t,x(:,4));
title('vitesse en z m/s');
xlabel('t (x)');
ylabel('vZ (m)');
grid on;


%% Trajectoire et plaque
% Trajectory and plate position at starting, middle, and end time.
%  From center of gravity, neglect thickness.  A plate of length $l$
%  rotate of angle $\theta$ around the center of gravity.  The plate
%  have two corners.
%   
% $$x_1=x_{cg}+0.5 l \cos(-\theta)$$
%
% $$x_2=x_{cg}-0.5 l \cos(-\theta)$$
%
% $$y_1=y_{cg}+0.5 l \sin(-\theta)$$
%
% $$y_2=y_{cg}-0.5 l \sin(-\theta)$$

if ( dim > 1 )
figure;
% start plate
xcoin(1)=x(1,1)+cosd(-x(1,5))*l/2.0;
xcoin(2)=x(1,1)-cosd(-x(1,5))*l/2.0;
ycoin(1)=x(1,3)+sind(-x(1,5))*l/2.0;
ycoin(2)=x(1,3)-sind(-x(1,5))*l/2.0;
% middle position plate
mid=round(length(xspan)/2);
xcoinm(1)=x(mid,1)+cosd(-x(mid,5))*l/2.0;
xcoinm(2)=x(mid,1)-cosd(-x(mid,5))*l/2.0;
ycoinm(1)=x(mid,3)+sind(-x(mid,5))*l/2.0;
ycoinm(2)=x(mid,3)-sind(-x(mid,5))*l/2.0;
% end position plate
fin=length(xspan);
xcoine(1)=x(fin,1)+cosd(-x(fin,5))*l/2.0;
xcoine(2)=x(fin,1)-cosd(-x(fin,5))*l/2.0;
ycoine(1)=x(fin,3)+sind(-x(fin,5))*l/2.0;
ycoine(2)=x(fin,3)-sind(-x(fin,5))*l/2.0;
plot(x(:,1),x(:,3),xcoin,ycoin,xcoinm,ycoinm,xcoine,ycoine);
title('position');
xlabel('X (m)');
ylabel('Z (m)');
daspect([1 1 1]);
grid on;
end

%% resultats sous forme de tableaux

taille=int16(tmax/deltat);
Z=zeros(taille,1);
X=zeros(taille,1);
U=zeros(taille,1);
W=zeros(taille,1);
OMEGA=zeros(taille,1);
THETA=zeros(taille,1);
for i = 1:(taille)
    X(i,1)=x(i,1);
    Z(i,1)=x(i,3);
    U(i,1)=x(i,2);
    W(i,1)=x(i,4);
    OMEGA(i,1)=x(i,6);
% keep theta between 0 and 180
    THETA(i,1)=mod(x(i,5),360);
%    if THETA(i,1) > 180
%        THETA(i,1)=360-THETA(i,1);
%    end
%    THETA(i,1)=x(i,5);
end


%%
%   exporter les tableaux dans excel
%filename = 'Magnus.xlsx';
%xlswrite(filename,X,4,'V2')
%xlswrite(filename,Z,4,'W2')
% xlswrite(filename,THETA,3,'P2')
% xlswrite(filename,U,3,'AG2')
% xlswrite(filename,W,3,'AI2')
% xlswrite(filename,OMEGA,3,'AK2')

%%
%   Modif 2 bis : Sortie txt pour compatibilite Linux
%   Precision a 10^-9
header1 = 'X';
header2 = 'Y';
header3 = 'U';
header4 = 'V';
header5 = 'Theta';
header6 = 'Omega';
T=table(X, Z, U, W, THETA, OMEGA);
writetable(T,'Magnus.txt','Delimiter',' ')
%fid=fopen('Magnus.txt','w');
%fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 ...
%    '\t' header5 '\t' header6 '\n']);
%fprintf(fid, '%11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \n', ...
%    [X Z U W THETA OMEGA]);
%fclose(fid);

end

