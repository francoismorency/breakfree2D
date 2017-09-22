function [CN_theta] = CN_2D(gamma)
% interpolate the 2D value for CN from data base in a file

% this should be done only once at calculation start
%read data file in to x
%the table angle range from 0 to 90
x=csvread('cn2d-alpha.csv',6,0);
xx=x(:,1);
v=x(:,2);
% build piece wise polynomial for interpolation
pp=griddedInterpolant(xx,v,'spline');
%pp = interp1(xx,v,'spline','pp')
%-------------------------------------------------


for i=1:length(gamma)
%transform gamma (0-360) between (0-90)
    alpha=mod(gamma(i),180);
    if alpha > 90
        alpha=180-alpha;
    end

%interpolation table from data
    CN_theta(i)=pp(alpha);
end
