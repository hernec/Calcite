function [Ee,Eo] = e_and_o_components(OA,normal,Einc)

%OA = optic axis, Ei = incident electric field, normal = normal to E-field
%plane

%calculate extraordinary axis
proj_e = OA - dot(OA,normal)*normal; %Project OA onto plane of face. E field to find extraordinary axis of face 1
%project projection on face onto xy plane, which is perp to plane of E
%field
xynorm = [0;0;1];
e_axis = proj_e - dot(proj_e,xynorm)*xynorm;
%project E onto projection of OA to find the extraordinary component
Ee = (dot(Einc,e_axis))*e_axis/(norm(e_axis))^2;

%calculate ordinary axis on crystal face as perpendicular to the plane formed by the optic
%axis and the extraordinary axis projection on crystal face
proj_o = cross(OA,proj_e);
%project projection on face onto xy plane
o_axis = proj_o - dot(proj_o,xynorm)*xynorm;
o_axis = 2*o_axis/norm(o_axis);

%project E onto projection of HI to find the extraordinary component
Eo = (dot(Einc,o_axis))*o_axis/(norm(o_axis))^2;

% plot3([0;e_axis(1)],[0;e_axis(2)],[0;e_axis(3)],'g','Linewidth',2)
% plot3([0;1*o_axis(1)],[0;1*o_axis(2)],[0;1*o_axis(3)],'b','Linewidth',2)
% plot3([0;Ee(1)],[0;Ee(2)],[0;Ee(3)],'r','Linewidth',2)
% plot3([0;Eo(1)],[0;Eo(2)],[0;Eo(3)],'b','Linewidth',2)
