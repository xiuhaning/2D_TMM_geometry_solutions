function [hex_coord,theta,alpha,gamma,phi]=solve_hexagon_analytical(theta0,alpha0,gamma0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global a_b; global b_b; global c_b;
global psi_ab; global psi_bb; global psi_cb;
global a_r; global b_r; global c_r;
global psi_ar; global psi_br; global psi_cr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A_x = 0;
A_y = 0;
F_x = b_r*cos(theta0+psi_bb+psi_cr);
F_y = -b_r*sin(theta0+psi_bb+psi_cr);
E_x = b_r*cos(theta0+psi_bb+psi_cr)-b_b*cos(theta0+alpha0+psi_bb+psi_cr+psi_cb+psi_ar); 
E_y = -b_r*sin(theta0+psi_bb+psi_cr)+b_b*sin(theta0+alpha0+psi_bb+psi_cr+psi_cb+psi_ar); 

B_x = c_b;
B_y = 0;
C_x = c_b-c_r*cos(gamma0+psi_ab+psi_br);
C_y = -c_r*sin(gamma0+psi_ab+psi_br);

A_angle = 2*pi-theta0-psi_bb-psi_cr; %Angle A found by subtracting known triangle angles
B_angle = 2*pi-gamma0-psi_ab-psi_br; %Angle F
F_angle = 2*pi-alpha0-psi_cb-psi_ar; %Angle B

%%
% Solving for D using analytical solutions
D_x1 = - (a_r^2 - E_x^2 - E_y^2 - a_b^2 + C_x^2 + C_y^2)/(2*(E_x - C_x)) - ((E_y - C_y)*(E_x^2*E_y - a_r^2*E_y + E_y*a_b^2 + a_r^2*C_y + E_y*C_x^2 + E_x^2*C_y - E_y*C_y^2 - E_y^2*C_y - a_b^2*C_y + C_x^2*C_y + E_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) - C_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + E_y^3 + C_y^3 - 2*E_x*E_y*C_x - 2*E_x*C_x*C_y))/(2*(E_x - C_x)*(E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y + C_x^2 + C_y^2));
D_x2 = - (a_r^2 - E_x^2 - E_y^2 - a_b^2 + C_x^2 + C_y^2)/(2*(E_x - C_x)) - ((E_y - C_y)*(E_x^2*E_y - a_r^2*E_y + E_y*a_b^2 + a_r^2*C_y + E_y*C_x^2 + E_x^2*C_y - E_y*C_y^2 - E_y^2*C_y - a_b^2*C_y + C_x^2*C_y - E_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + C_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + E_y^3 + C_y^3 - 2*E_x*E_y*C_x - 2*E_x*C_x*C_y))/(2*(E_x - C_x)*(E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y + C_x^2 + C_y^2));
D_y1 = (E_x^2*E_y - a_r^2*E_y + E_y*a_b^2 + a_r^2*C_y + E_y*C_x^2 + E_x^2*C_y - E_y*C_y^2 - E_y^2*C_y - a_b^2*C_y + C_x^2*C_y + E_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) - C_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + E_y^3 + C_y^3 - 2*E_x*E_y*C_x - 2*E_x*C_x*C_y)/(2*(E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y + C_x^2 + C_y^2));
D_y2 = (E_x^2*E_y - a_r^2*E_y + E_y*a_b^2 + a_r^2*C_y + E_y*C_x^2 + E_x^2*C_y - E_y*C_y^2 - E_y^2*C_y - a_b^2*C_y + C_x^2*C_y - E_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + C_x*((a_r^2 + 2*a_r*a_b - E_x^2 + 2*E_x*C_x - E_y^2 + 2*E_y*C_y + a_b^2 - C_x^2 - C_y^2)*(- a_r^2 + 2*a_r*a_b + E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y - a_b^2 + C_x^2 + C_y^2))^(1/2) + E_y^3 + C_y^3 - 2*E_x*E_y*C_x - 2*E_x*C_x*C_y)/(2*(E_x^2 - 2*E_x*C_x + E_y^2 - 2*E_y*C_y + C_x^2 + C_y^2));
%D1 is the concave solution, D2 is the convex solution

if imag(D_x1)~=0
   display('No solution for hexagon') %No D-point
   hex_coord = zeros(12,1);
   theta = 0;
   gamma = 0;
   alpha = 0;
   phi = 0;
   return 
end


%%
% solve for angles

%Choose either convex D1 or concave D2 
if theta0<pi
    D_x = D_x2;
    D_y = D_y2;
elseif theta0>pi
    D_x = D_x1;
    D_y = D_y1;
end

%solving for gamma using point D and E 
beta1=atan2(D_y-E_y,D_x-E_x);
if beta1<0
    beta1=beta1+2*pi;
end
beta2=2*pi-A_angle-F_angle;
if beta1+beta2<=2*pi
    gamma=beta1+beta2;
else
    gamma=beta1+beta2-2*pi;
end

%solving for alpha using point D and C 
delta_1=atan2(D_y-C_y,D_x-C_x);
if delta_1<0
    delta_1=delta_1+2*pi;
end
delta1=pi-delta_1;
delta2=pi-B_angle;
alpha=delta1+delta2;

theta=4*pi-A_angle-B_angle-F_angle-alpha-gamma;

if theta>2*pi-psi_cr-psi_bb || theta<0 || alpha>2*pi-psi_cb-psi_ar || alpha <0 || gamma> 2*pi-psi_ab-psi_br || gamma<0
   display('No solution for hexagon') %No D-point
   hex_coord = zeros(12,1);
   theta = 0;
   gamma = 0;
   alpha = 0;
   phi = 0;
   return 
end

phi=psi_bb-delta1;

hex_coord = [A_x A_y B_x B_y C_x C_y D_x D_y E_x E_y F_x F_y];


