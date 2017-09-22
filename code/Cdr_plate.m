function  [Cdr] = Cdr_plate(Urel,omega)

% corr�lation de Tachikawa valable pour un rapport l/e=0.05 pour la vitesse
% angulaire d'�quilibre
global l e L dim;
tau=e/l;
AR=L/l; %=largeur/longueur, donc 1 pour une plaque carree
So=(0.329*log(1/tau)-0.0246*(log(1/tau))^2)*((AR/(2+(4+AR^2)^0.5))...
    *(2-((AR)/(AR+0.595))^0.76))^(2/3);
% wo=0.64*Vf/l;

%Urel=((Vf-u)^2+(Wf-w)^2)^0.5; %S=omega*l/(2*Urel)
    
for i=1:length(omega)
S=omega(i)*l/(2*Urel);
if ( dim == 2)
    Cdro=1.28;
    if ( abs(S/So) <= 0.35 )
        Cdr = Cdro*(abs(S/So)+0.6);
    elseif ( 0.35 < abs(S/So) && abs(S/So) < 1 )
        Cdr = Cdro*(0.0769*(abs(S/So) - 0.35)) +0.95;
    else (abs(S/So) >= 1 )
        Cdr = Cdro;
    end
% mean value Cd statique 2D=1.1718
% see equation (5) Tachikawa(1983)
%"Trajectories of flat plates in uniform flow with 
%application to wind-generated missiles." 
%Journal of wind engineering and industrial aerodynamics 14: 443-453.
 
    Cdr = (Cdr - 1.1718/2.0);
else
    if (0.4 < (abs(S/So)) && (abs(S/So) < 1))
    
        Cdr=0.12+0.36*abs(S/So); %Cdr pour w/wo entre 0.4 et 1
    
    elseif abs(S/So)<= 0.4
    
        Cdr=0.66*abs(S/So); %Cdr pour w/wo inferieur a 0.4
    
    else
    
        Cdr=0.48; %Cdr pour w/wo superieur ou egal a 1

    end
end
end

% on ajoute ce coefficient au Cd deja existant pour Holmes.
% correlation de Journal of fluids and structures, D.M. Hargreaves,  B.
% Kakimpa, J.S. Owen 46 (2014) 111-133
