%Calcite_Torque.m
%C. Herne, 10/24/22

clear all

format short 


no = 1.658;
ne = 1.486;


nw = 1.33;


%Define constants
lambda = 660e-9; %units: m
c = 3e8; %m/s
omega= (2*pi*c/lambda);
e = (710e-12); %s^4/(kgm^3)A^2
e0 = 8.85e-12; %permittivity of free space, units F/m
k = (2*pi)/lambda;% 2pi/lambda  1/m?

d = 5*1e-6; % m,  max value for phi = pi/4 at d = 2.2 microns; 
%max for both spin and alignment is about d = 2.3 microns
Eo = 6e-3; % units?
P = 10e-3; %power; units of watts
C = 1; %This will be a constant of some magnitude, based on the above constants
R = 5e-6; %radius of crystal directions perpendicular to optic axis, meters
B = 12.6e-12; %overlap of beam with crystal in microns, for 2 micron radius get 12.6e-12

origin = [0;0;0];
E0 = 1;

%define incident ray (in z direction)
rayi = [0;0;1];
%define incident electric field
Ei = [E0;0;0];

%% Draw calcite and calculate face normals

alpha = 0;
%betaarray = [-(51+12)*pi/180:pi/36:(51-12)*pi/180];
betaarray = [0*pi/180:1*pi/180:1*pi/180];
gammaarray = [0*pi/180:1*pi/180:1*pi/180];

sb = size(betaarray);
Lb = sb(2);
sg = size(gammaarray);
Lg = sg(2);

% %plot calcite in zero position
% figure(1)
% D = 5; %in microns, for plotting
% [thetae1,thetae2,thetae3] = plotcalcite(D,0,0,0)


%hold on
% vectline(origin,P)
% vectline(origin,Ec)
% % plot torque
% vectline(origin,torque)

%hold off


torquevect1 = zeros();
torquevect2 = zeros();
torquevect3 = zeros();
ne_thetavect1 = zeros();
ne_thetavect2 = zeros();
ne_thetavect3 = zeros();
thetae1_array = zeros();
thetae2_array = zeros();
thetae3_array = zeros();

figure(1)
hold off

D = 5; %in microns, for plotting
[VertexArray] = plotcalcite_zero(D);
for i = 1:Lb
    for j = 1:Lg
    %plot calcite
    [thetae1,thetae2,thetae3,Ete1,Eto1,Ete2,Eto2,Ete3,Eto3,flux1,flux2,flux3,OA] = plotcalcite(VertexArray,alpha,betaarray(i),gammaarray(j),rayi,Ei);
    hold on
    [torque1,ne_theta1] = torque_PxE(thetae1,betaarray(i),gammaarray(j),Ete1,Eto1,ne,no,nw,1,flux1,OA);
    [torque2,ne_theta2] = torque_PxE(thetae2,betaarray(i),gammaarray(j),Ete2,Eto2,ne,no,nw,2,flux2,OA);
    [torque3,ne_theta3] = torque_PxE(thetae3,betaarray(i),gammaarray(j),Ete3,Eto3,ne,no,nw,3,flux3,OA);

    ne_thetavect1(i) = ne_theta1;
    ne_thetavect2(i) = ne_theta2;
    ne_thetavect3(i) = ne_theta3;
    thetae1_array(i) = thetae1;
    thetae2_array(i) = thetae2;
    thetae3_array(i) = thetae3;

    for k = 1:3
        torquevect1(i,j,k) = torque1(k);
        torquevect2(i,j,k) = torque2(k);
        torquevect3(i,j,k) = torque3(k);
        torquevecttotal(i,j,k) = torquevect1(i,j,k)+torquevect2(i,j,k)+torquevect3(i,j,k);
    end
    torquetotal(i,j) = norm(torquevecttotal(i,j,k));
    end
end

%% Find minimum torque value as beta changes for each gamma - find tilt angle
%vs. rotation angle for minimum torque.

for i = 1:Lg
    torquemin = torquetotal(1,i);
    true=0;
    false=0;
    count = i;
    betamin(i) = betaarray(1);
    betaminvalue = betamin(i);
    for j = 2:Lb
        if torquetotal(j,i) < torquemin
            torquemin = torquetotal(j,i);
            betaminvalue = betaarray(j);
            true = true+1;
        else
            false=false+1;
        end
    end
    torqueminarray(i) = torquemin;
    betamin(i) = betaminvalue;

end
save('betagammamin3.txt','betamin','gammaarray','torqueminarray','-ascii','-tabs')

%% Plotting

% figure(12)
% plot(gammaarray*180/pi,betamin*180/pi,'b*')
% xlabel('rotation about z axis (gamma)')
% ylabel('tilt (beta) for minimum torque')
% grid on
% g=gca;
% set(g,'fontsize',20)

figure(2)
%plot torque vs. incident angle - the last element of the array
plot(betaarray*180/pi,torquevect1(:,Lg-1,2),'r-.',betaarray*180/pi,torquevect2(:,Lg-1,2),'g--',betaarray*180/pi,torquevect3(:,Lg-1,2),'b:',betaarray*180/pi,torquevecttotal(:,Lg-1,2),'k','Linewidth',1.5)
%title('Y-direction torque on each face from PxE')
xlabel('rotation about y axis (beta)')
ylabel('torque in y-direction (arb. units)')
legend('section 3','section 2','section 1','sum')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)

figure(11)
%plot torque vs. incident angle - the last element of the array
plot(gammaarray*180/pi,torquevect1(Lb-1,:,2),'r-.',gammaarray*180/pi,torquevect2(Lb-1,:,2),'g--',gammaarray*180/pi,torquevect3(Lb-1,:,2),'b:',gammaarray*180/pi,torquevecttotal(Lb-1,:,2),'k','Linewidth',1.5)
%title('Y-direction torque on each face from PxE')
xlabel('rotation about z axis (gamma)')
ylabel('torque in y-direction (arb. units)')
legend('section 3','section 2','section 1','sum')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)

figure(3)
%plot total torque vs. incident angle
plot(betaarray*180/pi,torquevect1(:,Lg-1,1),'r-.',betaarray*180/pi,torquevect2(:,Lg-1,1),'g--',betaarray*180/pi,torquevect3(:,Lg-1,1),'b:',betaarray*180/pi,torquevecttotal(:,Lg-1,1),'k','Linewidth',1.5)
%title('X-direction torque on each face from PxE')
xlabel('rotation about y axis (beta)')
ylabel('torque in x-direction (arb. units)')
legend('face 1','face 2','face 3')
%xlim([0,4])
%ylim([0,1])
grid on
g=gca;
set(g,'fontsize',20)


figure(4)
%plot total torque vs. incident angle
plot(gammaarray*180/pi,torquevect1(Lb-1,:,3),'r-.',gammaarray*180/pi,torquevect2(Lb-1,:,3),'g--',gammaarray*180/pi,torquevect3(Lb-1,:,3),'b:',gammaarray*180/pi,torquevecttotal(Lb-1,:,3),'k','Linewidth',1.5)
%title('Z-direction torque on each face from PxE')
xlabel('rotation about z axis (gamma)')
ylabel('torque in z-direction (arb. units)')
legend('section 3','section 2','section 1','sum')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)

figure(5)
[X,Y] = meshgrid(gammaarray,betaarray);
%Cmap = [-1,0,1];
surf(X*180/pi,Y*180/pi,torquevecttotal(:,:,3))
colorbar
xlabel('gamma')
ylabel('beta')
title('total torque about z axis')

figure(6)
[X,Y] = meshgrid(gammaarray,betaarray);
surf(X*180/pi,Y*180/pi,torquevecttotal(:,:,2))
colorbar
xlabel('gamma')
ylabel('beta')
title('total torque about y axis')

figure(7)
[X,Y] = meshgrid(gammaarray,betaarray);
surf(X*180/pi,Y*180/pi,torquevecttotal(:,:,1))
colorbar
xlabel('gamma')
ylabel('beta')
title('total torque about x axis')

figure(8)
[X,Y] = meshgrid(gammaarray,betaarray);
surf(X*180/pi,Y*180/pi,torquetotal(:,:))
colorbar
xlabel('gamma')
ylabel('beta')
title('total torque norm for all axes')

figure(9)
plot(betaarray*180/pi,torquetotal(:,Lg-1),'k',betaarray*180/pi,torquevecttotal(:,Lg-1,1),betaarray*180/pi,torquevecttotal(:,Lg-1,2),betaarray*180/pi,torquevecttotal(:,Lg-1,3),'Linewidth',1.5)
%title('total torque on all faces')
xlabel('rotation about y axis (beta)')
ylabel('total torque on all faces (arb. units)')
legend('sum','x','y','z')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)

figure(12)
plot(betaarray*180/pi,torquetotal(:,Lg-1),'k','Linewidth',1.5)
%title('total torque (arb. units)')
xlabel('rotation about y axis (degrees)')
ylabel('total torque (arb. units)')
%legend('sum','x','y','z')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)

figure(10)
plot(gammaarray*180/pi,torquetotal(Lb-1,:),'k','Linewidth',1.5)
%title('total torque on all faces')
xlabel('rotation about z axis (degrees)')
ylabel('total torque (arb. units)')
%legend('sum')
%xlim([0,4])
%ylim([0,4])
grid on
g=gca;
set(g,'fontsize',20)
