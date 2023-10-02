function [torque,ne_theta] = torque_PxE(thetae,beta,gamma,Ete,Eto,ne,no,nw,face,flux,OA)

ne_theta = sqrt(1/(cos(thetae)^2/no^2 + sin(thetae)^2/ne^2)); %angle between propagation and optic axis
%Ete and Eto are transmitted e and o rays; transmitted ray in lab frame
% Ete
% Eto
Et = (Ete+Eto);
Etodir = Eto/norm(Eto);
Etedir = Ete/norm(Ete);
OAdir = OA/norm(OA);
%Define direction of polarization vector in the crystal. Switch y direction
%for -x direction because of how crystal is drawn in plotcalcite

%try defining directions of Pe based on Ete
P_e_x = (sin(thetae)/(ne_theta^2-no^2))*Etedir;
P_e_z = cos(thetae)/(ne_theta^2-ne^2)*OAdir;
P_e_temp = P_e_x+P_e_z;
% P_e_temp = [sin(thetae)/(ne_theta^2-no^2);0;cos(thetae)/(ne_theta^2-ne^2)]; %in crystal frame
P_e = P_e_temp/norm(P_e_temp);
P_o = Etodir;
P_o_mag = (no^2-1)*P_o*norm(Eto); %chi = n^2 - 1; P = chi E; multiply by magnitude of E field
P_e_mag = (ne^2-1)*P_e*norm(Ete);

%normalize
% Punit = P/norm(P);
% Pval = [(no^2-1)*Punit(1);(ne_theta^2-1)*Punit(2);(ne_theta^2-1)*Punit(3)];

%Ryneg = [cos(-beta) 0 sin(-beta); 0 1 0; -sin(-beta) 0 cos(-beta)];

if face == 1
 %       gamma1 = gamma - (180-53.7)*pi/180; %-(180-53.7)*pi/180; 2.38 rad
%     P_o = [0;-1;0];
%     P_o_mag = (no^2-1)*P_o*norm(Eto); %chi = n^2 - 1; P = chi E; multiply by magnitude of E field
%     P_e_mag = (ne^2-1)*P_e*norm(Ete); %chi = n^2 - 1; P = chi E
%     dP = abs(atan(norm(P_e_mag)/norm(P_o_mag)));
%     dPtest = dP*180/pi
%     dgamma = 103*pi/180; %103
%     gammatemp = gamma-dgamma;
%     gammatemptest = gammatemp*180/pi
% 
%     if -5*pi/2<=gammatemp & gammatemp<-2*pi
%         gamma1 = gammatemp - 2*(pi/2-dP) + 2*pi;
%     elseif -2*pi<=gammatemp & gammatemp<-3*pi/2
%         gamma1 = gammatemp - pi + pi;
%     elseif -3*pi/2<=gammatemp & gammatemp<-pi
%         gamma1 = gammatemp - 2*(pi/2-dP) + pi;
%     elseif -pi<=gammatemp & gammatemp<-pi/2
%         gamma1 = gammatemp - pi
%     elseif -pi/2<=gammatemp & gammatemp<0
%         gamma1 = gammatemp - 2*(pi/2-dP)
%     elseif 0<=gammatemp & gammatemp<pi/2
%         gamma1 = gammatemp;
%     elseif pi/2<=gammatemp & gammatemp<pi
%         gamma1 = gammatemp+2*dP;
%     elseif pi<=gammatemp & gammatemp<3*pi/2
%         gamma1 = gammatemp - pi;
%     elseif 3*pi/2<=gammatemp & gammatemp<=2*pi
%         gamma1 = gammatemp + 2*dP - pi;
%     elseif 2*pi<=gamma & gamma<5*pi/2
%         gamma1 = gammatemp - 2*pi;
%     elseif 5*pi/2<=gamma & gamma<=3*pi
%         gamma1 = gammatemp + 2*dP - 2*pi;
%     else gamma1 = gammatemp;
%     end
%     Rz = [cos(gamma1) -sin(gamma1) 0; sin(gamma1) cos(gamma1) 0; 0 0 1];
%     
%     if -pi/2<=beta & beta<0
%         beta1 = 1*beta;
%     elseif 0<=beta & beta<pi/2
%         beta1 = beta;
%     end
%     Ry = [cos(beta1) 0 sin(beta1); 0 1 0; -sin(beta1) 0 cos(beta1)];

    PE = P_e_mag + P_o_mag;
%     PErot = Rz*Ry*PE; %rotate PE into lab frame from 
    torque = cross(PE,Et)*flux;
    
% %      figure(1)
%             plot3([0;5*Et(1)],[0;5*Et(2)],[0;5*Et(3)],'c','Linewidth',1.5)
%       plot3([0;5*Ete(1)],[0;5*Ete(2)],[0;5*Ete(3)],'r','Linewidth',2)
%       plot3([0;5*Eto(1)],[0;5*Eto(2)],[0;5*Eto(3)],'m','Linewidth',2)
% % 
% %     plot3([0;PE(1)],[0;PE(2)],[0;PE(3)],'k','Linewidth',1.5)
% % %     plot3([0;PErot(1)],[0;PErot(2)],[0;PErot(3)],'r','Linewidth',1.5)
% %         plot3([0;1*P_e_mag(1)],[0;1*P_e_mag(2)],[0;1*P_e_mag(3)],'c','Linewidth',3)
% %     plot3([0;1*P_o_mag(1)],[0;1*P_o_mag(2)],[0;1*P_o_mag(3)],'g','Linewidth',3)
% % % 
% %     plot3([0;10*torque(1)],[0;10*torque(2)],[0;10*torque(3)],'m','Linewidth',1.5)  

elseif face == 2
    
%     P_o = [0;-1;0];
%     P_o_mag = (no^2-1)*P_o*norm(Eto); %chi = n^2 - 1; P = chi E; multiply by magnitude of E field
% %     P_e_mag = (ne^2-1)*P_e*norm(Ete); %chi = n^2 - 1; P = chi E
%     dP = abs(atan(norm(P_e_mag)/norm(P_o_mag)));
%     dgamma = 250*pi/180;
%     gammatemp = gamma-dgamma;
% 
%     if -5*pi/2<=gammatemp & gammatemp<-2*pi
%         gamma2 = gammatemp - 2*(pi/2-dP) + 2*pi;
%     elseif -2*pi<=gammatemp & gammatemp<-3*pi/2
%         gamma2 = gammatemp - pi + pi;
%     elseif -3*pi/2<=gammatemp & gammatemp<-pi
%         gamma2 = gammatemp - 2*(pi/2-dP) + pi;
%     elseif -pi<=gammatemp & gammatemp<-pi/2
%         gamma2 = gammatemp - pi;
%     elseif -pi/2<=gammatemp & gammatemp<0
%         gamma2 = gammatemp - 2*(pi/2-dP);
%     elseif 0<=gammatemp & gammatemp<pi/2
%         gamma2 = gammatemp;
%     elseif pi/2<=gammatemp & gammatemp<pi
%         gamma2 = gammatemp+2*dP;
%     elseif pi<=gammatemp & gammatemp<3*pi/2
%         gamma2 = gammatemp - pi;
%     elseif 3*pi/2<=gammatemp & gammatemp<=2*pi
%         gamma2 = gammatemp + 2*dP - pi;
%     elseif 2*pi<=gamma & gamma<5*pi/2
%         gamma2 = gammatemp - 2*pi;
%     elseif 5*pi/2<=gamma & gamma<=3*pi
%         gamma2 = gammatemp + 2*dP - 2*pi;
%     else gamma2 = gammatemp;
%     end
%     Rz = [cos(gamma2) -sin(gamma2) 0; sin(gamma2) cos(gamma2) 0; 0 0 1];
%     
%     if -pi/2<=beta & beta<0
%         beta2 = beta;
%     elseif 0<=beta & beta<pi/2
%         beta2 = beta;
%     end
%      Ry = [cos(beta2) 0 sin(beta2); 0 1 0; -sin(beta2) 0 cos(beta2)];

    PE = P_e_mag + P_o_mag;
%     PErot = Rz*Ry*PE;
    torque = cross(PE,Et)*flux;

elseif face == 3
    
%  
%     dP = abs(atan(norm(P_e_mag)/norm(P_o_mag)));
%     dPy = abs(atan((P_e_mag(3))/(P_e_mag(1))));
%    
%     if -2*pi<=gamma & gamma<-3*pi/2
%         gamma3 = gamma - pi + pi;
%     elseif -3*pi/2<=gamma & gamma<-pi
%         gamma3 = gamma - 2*(pi/2-dP) + pi;
%     elseif -pi<=gamma & gamma<-pi/2
%         gamma3 = gamma - pi;
%     elseif -pi/2<=gamma & gamma<0
%         gamma3 = gamma - 2*(pi/2-dP);
%     elseif 0<=gamma & gamma<pi/2
%         gamma3 = gamma;
%     elseif pi/2<=gamma & gamma<pi
%         gamma3 = gamma+2*dP;
%     elseif pi<=gamma & gamma<3*pi/2
%         gamma3 = gamma - pi;
%     elseif 3*pi/2<=gamma & gamma<=2*pi
%         gamma3 = gamma + 2*dP - pi;
%     elseif 2*pi<=gamma & gamma<5*pi/2
%         gamma3 = gamma - 2*pi;
%     elseif 5*pi/2<=gamma & gamma<=3*pi
%         gamma3 = gamma + 2*dP - 2*pi;
%     end
%      Rz = [cos(gamma3) -sin(gamma3) 0; sin(gamma3) cos(gamma3) 0; 0 0 1]; 
%      
%      if -pi/2<=beta & beta<0
%         beta3 = beta;
%     elseif 0<=beta & beta<pi/2
%         beta3 = beta;
%      end
%      Ry = [cos(beta3) 0 sin(beta3); 0 1 0; -sin(beta3) 0 cos(beta3)];

    PE = P_e_mag + P_o_mag;
%     PErot = Rz*Ry*PE;
    torque = cross(PE,Et)*flux;

plot3([0;4*PE(1)],[0;4*PE(2)],[0;4*PE(3)],'g','Linewidth',1.5)
plot3([0;3*Et(1)],[0;3*Et(2)],[0;3*Et(3)],'r','Linewidth',1.5)


end

