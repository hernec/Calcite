%Draw_calcite.m


clear all


no = 1.658;
ne = 1.486;
nw = 1.33;

origin = [0;0;0];
E0 = 1;

%define incident ray (in z direction)
rayi = [0;0;1];
%define incident electric field
Ei = [E0;0;0];

D = 5; %in microns, for plotting

%% Draw zero position

[vertices] = plotcalcite_zero(D);

Atemp = vertices(:,1);
Btemp = vertices(:,2);
Ctemp = vertices(:,3);
Dtemp = vertices(:,4);
Etemp = vertices(:,5);
Ftemp = vertices(:,6);
Gtemp = vertices(:,7);
Htemp = vertices(:,8);
Itemp = vertices(:,9);

%% now rotate by angles 
alpha0 = 0*pi/180;
beta = 0*pi/180;
gamma = 20*pi/180;

%Enter rotation matrices 
    Rx = [1 0 0; 0 cos(alpha0) -sin(alpha0); 0 sin(alpha0) cos(alpha0)];
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


%calculate extraordinary axis for each face
%make vector that is along optic axis
HI = (Irot-Hrot);
OA = HI;

figure(1)
hold on
axis equal
%calculate extraordinary components
[Ee1,Eo1] = e_and_o_components(HI,normal1,Ei);
[Ee2,Eo2] = e_and_o_components(HI,normal2,Ei);
[Ee3,Eo3] = e_and_o_components(HI,normal3,Ei);

%Draw lines around the entire rhombohedron shape
    xarray = [Arot(1) Brot(1) Crot(1) Drot(1) Hrot(1) Grot(1) Crot(1) Brot(1) Frot(1) Erot(1) Arot(1) Drot(1) Hrot(1) Erot(1) Frot(1) Grot(1)];
    yarray = [Arot(2) Brot(2) Crot(2) Drot(2) Hrot(2) Grot(2) Crot(2) Brot(2) Frot(2) Erot(2) Arot(2) Drot(2) Hrot(2) Erot(2) Frot(2) Grot(2)];
    zarray = [Arot(3) Brot(3) Crot(3) Drot(3) Hrot(3) Grot(3) Crot(3) Brot(3) Frot(3) Erot(3) Arot(3) Drot(3) Hrot(3) Erot(3) Frot(3) Grot(3)];

%% Plot calcite

figure(1)
plot3(xarray, yarray, zarray,'b-','Linewidth',1.5)
hold on
xlabel('x')
ylabel('y')
zlabel('z')
%grid on
axis ([-5,5, -5, 5, -5 ,5]) 
plot3([Hrot(1);Irot(1)],[Hrot(2);Irot(2)],[Hrot(3);Irot(3)],'r','Linewidth',1.5)
plot3([-5;5],[0;0],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[-5;5],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[0;0],[-5;5],'k','Linewidth',1.5') %Make axes lines
% plot3([0;thetat3vect(1)],[0;thetat3vect(2)],[0;thetat3vect(3)],'r')
% plot3([0;thetat2vect(1)],[0;thetat2vect(2)],[0;thetat2vect(3)],'g')
% plot3([0;thetat1vect(1)],[0;thetat1vect(2)],[0;thetat1vect(3)],'b')
%  plot3([0;1*Eto1(1)],[0;1*Eto1(2)],[0;1*Eto1(3)],'r','Linewidth',2)
%  plot3([0;Ete1(1)],[0;Ete1(2)],[0;Ete1(3)],'m','Linewidth',2)
set(gca,'fontsize',15)
axis equal

%% Plot faces

X1 = [Arot(1) Brot(1) Crot(1) Drot(1)];
Y1 = [Arot(2) Brot(2) Crot(2) Drot(2)];
Z1 = [Arot(3) Brot(3) Crot(3) Drot(3)];

X2 = [Brot(1) Crot(1) Grot(1) Frot(1)];
Y2 = [Brot(2) Crot(2) Grot(2) Frot(2)];
Z2 = [Brot(3) Crot(3) Grot(3) Frot(3)];

X3 = [Crot(1) Drot(1) Hrot(1) Grot(1)];
Y3 = [Crot(2) Drot(2) Hrot(2) Grot(2)];
Z3 = [Crot(3) Drot(3) Hrot(3) Grot(3)];

X4 = [Erot(1) Frot(1) Grot(1) Hrot(1)];
Y4 = [Erot(2) Frot(2) Grot(2) Hrot(2)];
Z4 = [Erot(3) Frot(3) Grot(3) Hrot(3)];

X5 = [Arot(1) Drot(1) Hrot(1) Erot(1)];
Y5 = [Arot(2) Drot(2) Hrot(2) Erot(2)];
Z5 = [Arot(3) Drot(3) Hrot(3) Erot(3)];

X6 = [Arot(1) Brot(1) Frot(1) Erot(1)];
Y6 = [Arot(2) Brot(2) Frot(2) Erot(2)];
Z6 = [Arot(3) Brot(3) Frot(3) Erot(3)];

XP1 = [Brot(1) Frot(1) Hrot(1) Drot(1)];
YP1 = [Brot(2) Frot(2) Hrot(2) Drot(2)];
ZP1 = [Brot(3) Frot(3) Hrot(3) Drot(3)];

%top faces
fill3(X1,Y1,Z1,'b'); alpha .3
fill3(X2,Y2,Z2,'b'); alpha .3
fill3(X6,Y6,Z6,'b'); alpha .6
%bottom faces
fill3(X3,Y3,Z3,'b'); alpha 0.3
fill3(X4,Y4,Z4,'b'); alpha 0.3
fill3(X5,Y5,Z5,'b'); alpha 0.3

%principal section
fill3(XP1,YP1,ZP1,'r'); alpha 0.3

hold off
%% Plot sphere

 %plot ellipsoid
 [X,Y,Z] = ellipsoid(0,0,0,1,1,1);
 figure(2)
 hold on
 surf(X,Y,Z)
 colormap('default')
 plot3([-1.5;1.5],[0;0],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[-1.5;1.5],[0;0],'k','Linewidth',1.5) %Make axes lines
plot3([0;0],[0;0],[-1.5;1.5],'k','Linewidth',1.5') %Make axes lines
 %colorbar
 xlabel('x')
 ylabel('y')
 zlabel('z')
 axis equal
 ax = gca;
 alpha 0.5
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'fontsize',20)
%set(gca,'XTick',[], 'YTick', [],'ZTick',[]) %hides axes values