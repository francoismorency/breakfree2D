function [alpha] = alpha(u, w, vf, wf)
% Compute the relative velocity of the particule

% first, compute relative velocity component
urel=vf-u;
% this sign of the w freestream velocity must be coherent with gravity 
wrel=wf-w;

% angle du vent relatif en degr�s

alpha =atan2(wrel,urel)*180.0/pi;

end
