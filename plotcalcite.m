function [thetae1e,thetae2e,thetae3e,Ete1,Eto1,Ete2,Eto2,Ete3,Eto3,flux1,flux2,flux3,OA] = plotcalcite(vertices,alpha,beta,gamma,rayi,Ei)
%plots a calcite rhomb with edge length d
%define variables
%alpha is rotation about x, beta is rotation about y, gamma is rotation
%about z
% 
% no = 5;
% ne = 1.5;
no = 1.658;
ne = 1.486;
nw = 1.33;

k=1;%value of proportion of force to intensity
lambda=0.6660*1e0; %wavelength
m=1; %maximum intensity
w=1e0; %minimum size of beam waist

%%
Atemp = vertices(:,1);
Btemp = vertices(:,2);
Ctemp = vertices(:,3);
Dtemp = vertices(:,4);
Etemp = vertices(:,5);
Ftemp = vertices(:,6);
Gtemp = vertices(:,7);
Htemp = vertices(:,8);
Itemp = vertices(:,9);

%now rotate by angle given in input arguments

%Enter rotation matrices 
    Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
    Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
    Ryneg = [cos(-beta) 0 sin(-beta); 0 1 0; -sin(-beta) 0 cos(-beta)];
    Rz = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
    Rzneg = [cos(-gamma) -sin(-gamma) 0; sin(-gamma) cos(-gamma) 0; 0 0 1];

    %rotate A through I vectors by angles, first about z then x then y
    Arot = Rx*Rz*Ry*Atemp;
    Brot = Rx*Rz*Ry*Btemp;
    Crot = Rx*Rz*Ry*Ctemp;
    Drot = Rx*Rz*Ry*Dtemp;
    Erot = Rx*Rz*Ry*Etemp;
    Frot = Rx*Rz*Ry*Ftemp;
    Grot = Rx*Rz*Ry*Gtemp;
    Hrot = Rx*Rz*Ry*Htemp;
    Irot = Rx*Rz*Ry*Itemp;
    
%Calculate face normals and then normalize the vector and point it into the
%crystal
%face1 = HDAE, face2 = HDCG, face3 = HGFE

vec1 = cross(Hrot-Arot,Drot-Erot);
vec2 = cross(Hrot-Crot,Grot-Drot);
vec3 = cross(Hrot-Frot,Erot-Grot);
normal1 = vec1/norm(vec1);
normal2 = vec2/norm(vec2);
normal3 = vec3/norm(vec3);

%calculate flux through each face, using flux = A cos(theta) where A is
%area of face. All faces have same area so let A = 1. Then theta can be
%found by norm x B = norm B sin(theta)
planeperp = [0;0;1];
flux1 = abs(dot(normal1,planeperp));
flux2 = abs(dot(normal2,planeperp));
flux3 = abs(dot(normal3,planeperp));

% flux1 = cos(pi/2 - asin(norm(cross(normal1,plane))))
% flux2 = cos(pi/2 - asin(norm(cross(normal2,plane))))
% flux3 = cos(pi/2 - asin(norm(cross(normal3,plane))))


%calculate extraordinary axis for each face
%make vector that is along optic axis
HI = (Irot-Hrot);
OA = HI;

[Ee1,Eo1] = e_and_o_components(HI,normal1,Ei);
[Ee2,Eo2] = e_and_o_components(HI,normal2,Ei);
[Ee3,Eo3] = e_and_o_components(HI,normal3,Ei);

%calculate transmitted vector in crystal
[thetat1e,thetat1vect,Ete1] = thetatcalc(Ee1,rayi,normal1,gamma,nw,ne,1);
[thetat2e,thetat2vect,Ete2] = thetatcalc(Ee2,rayi,normal2,gamma,nw,ne,2);
[thetat3e,thetat3vect,Ete3] = thetatcalc(Ee3,rayi,normal3,gamma,nw,ne,3);
[thetat1o,thetat1vecto,Eto1] = thetatcalc(Eo1,rayi,normal1,gamma,nw,no,1);
[thetat2o,thetat2vecto,Eto2] = thetatcalc(Eo2,rayi,normal2,gamma,nw,no,2);
[thetat3o,thetat3vecto,Eto3] = thetatcalc(Eo3,rayi,normal3,gamma,nw,no,3);
Et1 = Ete1+Eto1;
Et2 = Ete2+Eto2;
Et3 = Ete3+Eto3;
%calculate angle between transmitted vector and optic axis
thetae1e = asin(norm(cross(thetat1vect,HI))/(norm(thetat1vect)*norm(HI)));
thetae2e = asin(norm(cross(thetat2vect,HI))/(norm(thetat2vect)*norm(HI)));
thetae3e = asin(norm(cross(thetat3vect,HI))/(norm(thetat3vect)*norm(HI)));
%thetae1test = thetae1e*180/pi

    
%Draw lines around the entire rhombohedron shape
    xarray = [Arot(1) Brot(1) Crot(1) Drot(1) Hrot(1) Grot(1) Crot(1) Brot(1) Frot(1) Erot(1) Arot(1) Drot(1) Hrot(1) Erot(1) Frot(1) Grot(1)];
    yarray = [Arot(2) Brot(2) Crot(2) Drot(2) Hrot(2) Grot(2) Crot(2) Brot(2) Frot(2) Erot(2) Arot(2) Drot(2) Hrot(2) Erot(2) Frot(2) Grot(2)];
    zarray = [Arot(3) Brot(3) Crot(3) Drot(3) Hrot(3) Grot(3) Crot(3) Brot(3) Frot(3) Erot(3) Arot(3) Drot(3) Hrot(3) Erot(3) Frot(3) Grot(3)];


plot3(xarray, yarray, zarray,'b-','Linewidth',1.5)
hold on
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis ([-4,4, -4, 4, -4 ,4]) 
plot3([Hrot(1);Irot(1)],[Hrot(2);Irot(2)],[Hrot(3);Irot(3)],'r')
plot3([-5;5],[0;0],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[-5;5],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[0;0],[-5;5],'k','Linewidth',1.5') %Make axes lines
% plot3([0;thetat3vect(1)],[0;thetat3vect(2)],[0;thetat3vect(3)],'r')
% plot3([0;thetat2vect(1)],[0;thetat2vect(2)],[0;thetat2vect(3)],'g')
% plot3([0;thetat1vect(1)],[0;thetat1vect(2)],[0;thetat1vect(3)],'b')
%  plot3([0;1*Eto1(1)],[0;1*Eto1(2)],[0;1*Eto1(3)],'r','Linewidth',2)
%  plot3([0;Ete1(1)],[0;Ete1(2)],[0;Ete1(3)],'m','Linewidth',2)

set(gca,'fontsize',20)


axis equal
%plot3([0;thetat3vect(1)],[0;thetat3vect(2)],[0;thetat3vect(3)],'k')
hold off