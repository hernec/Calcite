function [thetat,thetatvect,Et] = thetatcalc(Ei,incident,normal,gamma,nw,n,face)
%find transmitted ray and transmitted E field vector


thetai = asin(norm(cross(incident,normal))); %angle between incident ray and face normal

thetat = asin(nw*sin(thetai)/n); %angle between face normal and transmitted ray

perptemp = cross(incident,normal);  %perpendicular to plane containing incident and normal rays
perp = perptemp/norm(perptemp);
%rotate incident vector by thetat
A = [normal(1),normal(2),normal(3);perp(1),perp(2),perp(3);incident(1),incident(2),incident(3)];
B = [cos(thetat);0;cos(thetai-thetat)];
X = linsolve(A,B); %solve system of linear equations

thetatvect = X;
%%
%now rotate Ei so it's the transmitted ray Et
%rotate by thetai - thetat
% if face == 1
%     if normal(1) < 0 %& normal(1) < 0
%     theta = (thetai-thetat);
%     else theta = -(thetai-thetat);
%     end
% elseif face == 2
%     if normal(2) < 0 | normal(1) < 0
%          theta = -(thetai-thetat);
%      else theta = (thetai-thetat);
%      end
% elseif face == 3
%     if normal(2) < 0  | normal(1) < 0
%          theta = -(thetai-thetat);
%      else theta = (thetai-thetat);
%     end
% end 

theta = -(thetai-thetat);

%define unit vector (a,b,c) along the rotation axis
a = perp(1); b = perp(2); c = perp(3);
d = sqrt(b^2+c^2);
%define rotation matrices
Rx = [1,0,0,0;0,c/d,-b/d,0;0,b/d,c/d,0;0,0,0,1];
Rxinv = [1,0,0,0;0,c/d,b/d,0;0,-b/d,c/d,0;0,0,0,1];
Ry = [d,0,-a,0;0,1,0,0;a,0,d,0;0,0,0,1];
Ryinv = [d,0,a,0;0,1,0,0;-a,0,d,0;0,0,0,1];
Rz = [cos(theta),sin(theta),0,0;-sin(theta),cos(theta),0,0;0,0,1,0;0,0,0,1];


%operate on Ei

Et_temp = Rxinv*Ryinv*Rz*Ry*Rx*[Ei(1);Ei(2);Ei(3);1];
Et_temp_3D = [Et_temp(1);Et_temp(2);Et_temp(3)];
% Rzgamma = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
% 
% Et = Rzgamma*Et_temp_3D;
Et = Et_temp_3D;

%check = dot(Et,thetatvect)
