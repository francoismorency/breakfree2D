function  [Clr] = Clr_plate(Urel,omega)

% Tachikawa equations
% Calculation of the additional lift coefficient due to the Magnus effect
% in term of the angular velocity.
%Kordi, B. and G. A. Kopp (2009). "Evaluation of quasi-steady theory 
%applied to windborne flat plates in uniform flow." 
%Journal of engineering mechanics 135(7).
% 
global l e L dim;
tau=e/l;
AR=L/l; %=largeur/longueur, donc 1 pour une plaque carr�e
So=(0.329*log(1/tau)-0.0246*(log(1/tau))^2)*((AR/(2+(4+AR^2)^0.5))...
    *(2-((AR)/(AR+0.595))^0.76))^(2/3);
% wo=0.64*Vf/l;

%Urel=((Vf-u)^2+(Wf-w)^2)^0.5; %S=omega*l/(2*Urel)
    
 for i=1:length(omega)
    S=omega(i)*l/(2*Urel);
    if (dim == 2)
        Clro=0.76;
        if ( S/So <= -0.6 )
            Clr = Clro*((S/So+1)-1.0);
        elseif ( S/So > -0.6 && S/So <= -0.1 )
            Clr = Clro*((S/So+0.6)*0.5-0.6);
        elseif ( S/So > -0.1 && S/So <= -0.1 )
            Clr = Clro*(3.5*S/So);
        elseif ( S/So > -0.1 && S/So <= 0.6)
            Clr = Clro*(0.5*(S/So-0.1)+0.1);
        elseif (S/So > 0.6 )
            Clr = Clro*((S/So-0.6)+0.6);
        end
    else    
            
        if  (-0.2 < S /So) && (S/So < 0.2)
    
            Clr=1.05*(S/So); %Clr pour w/wo entre -0.2 et 0.2
    
        elseif ((S/So)<=-0.2)
    
            Clr=-0.1575+0.2625*(S/So); %Clr pour w/wo inferieur a -0.2
    
        else
    
            Clr=0.1575+0.2625*(S/So); %Clr pour w/wo superieur ou egal a 0.2

        end
    end
end


% on ajoute ce coefficient au Cl deja existant pour Holmes.
%correlation de Journal of fluids and structures, D.M. Hargreaves,  B.
%Kakimpa, J.S. Owen 46 (2014) 111-133



