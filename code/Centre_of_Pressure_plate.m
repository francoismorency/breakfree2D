function  [Lx] = Centre_of_Pressure_plate(gamma)

% Holmes equations
% Calcul de la positon du centre de pression par rapport au centre de la
% plaque
% Lx=c/l 
% c=distance du centre de pression au centre de la plaque
% l=corde de la plaque 
% gamma angle d'attaque total

gammad=mod(gamma,180);
for i=1:length(gamma)

  if (gammad(i) <= 38)                               %%When theta<=38�            
    Lx(i) = 0.3-0.22*(gammad(i)./38)    ;            %%Eq.42a
    
  elseif (38 < gammad(i)) && (gammad(i) < 87)        %%When 38<theta<82.5� 
    Lx(i) = 0.08.*cosd(2.*(gammad(i)-38))  ;        %%Eq.42b
    
  elseif (87 <= gammad(i)) && (gammad(i) < 92)    %%When 82.5<=theta<=97.5� 
    Lx(i) = 0;                                      %%Eq.42c
    
  elseif (92 <= gammad(i)) && (gammad(i) <= 142)    %%When 97.5<=theta<=142� 
    Lx(i) = -0.08*cosd(2*(142-gammad(i)))        ;   %%Eq.42d
    
  elseif (142 < gammad(i)) && (gammad(i) <= 181)     %%When 142<theta<=180 
    Lx(i) = -0.3+0.22*((180-gammad(i))./38)      ;   %%Eq.42e
    
       
  end
  

end
end