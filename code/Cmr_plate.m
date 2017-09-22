function  [Cmr] = Cmr_plate(Urel,omega,u,w)
%%
% Magnus effect for plate
global l L e dim
%l=0.04;
%e=0.002;
tau=e/l;
AR=L/l; %=largeur/longueur, donc 1 pour une plaque carr�e
So=(0.329*log(1/tau)-0.0246*(log(1/tau))^2)*((AR/(2+(4+AR^2)^0.5))*(2-((AR)/(AR+0.595))^0.76))^(2/3);
% wo=0.64*Vf/l;

%Urel=((Vf-u)^2+(Wf-w)^2)^0.5; %S=omega*l/(2*Urel)

for i=1:length(omega)
S= (omega(i)*l/(2*Urel)) ; 
if (dim == 2)
    if ( abs(S/So) <= 0.40 )
        Cmr = sign(S)*0.3*abs(S/So)-0.04;
    elseif ( 0.40 < abs(S/So) <= 0.70 )
        Cmr = sign(S)*0.08;
    elseif ( abs(S/So) > 0.70 )
        Cmr = sign(S)*(-0.2667*(abs(S/So)-0.7)+0.08);
    end
else
    
    if S/So <= 1 && S/So >= -1
    
    Cmr=0.12*(1-(abs(S/So)))*(S/So); %Cdr pour w/wo entre -1 et 1
    
    elseif S/So < -1
    
    Cmr=-0.12*(1+(S/So)); %Cmr pour w/wo inferieur a -1
    
    else
    
    Cmr=0.12*(1-(S/So)); %Cmr pour w/wo superieur  a 1

    end
end
end

% on ajoute ce coefficient au Cm deja existant pour Holmes.
% correlation de Journal of fluids and structures, D.M. Hargreaves,  B.
% Kakimpa, J.S. Owen 46 (2014) 111-133, or Tachikawa
