%% Solves Maxwell Lattice Soft Edge


%% Scrub-a-dub
clear all
clc
% close all
%% load Homogeneous lattice angles
load Homogeneous_lattice_angles.mat

%% Known Values
%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
n=5;m=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blue Triangle
%Define dimensions of blue triangle
global a_b; global b_b; global c_b;
%a_b=0.5;b_b=0.7;c_b=1;
a_b=0.8;b_b=1;c_b=0.4;

global psi_ab; global psi_bb; global psi_cb;
psi_ab=acos((a_b^2-b_b^2-c_b^2)/(-2*b_b*c_b));
psi_bb=acos((b_b^2-a_b^2-c_b^2)/(-2*a_b*c_b));
psi_cb=acos((c_b^2-b_b^2-a_b^2)/(-2*b_b*a_b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Red Triangle
%Define dimensions of red triangle
global a_r; global b_r; global c_r;
%a_r=0.4;b_r=0.8;c_r=1;
a_r=1;b_r=0.5;c_r=0.7;

global psi_ar; global psi_br; global psi_cr;
psi_ar=acos((a_r^2-b_r^2-c_r^2)/(-2*b_r*c_r));
psi_br=acos((b_r^2-a_r^2-c_r^2)/(-2*a_r*c_r));
psi_cr=acos((c_r^2-b_r^2-a_r^2)/(-2*b_r*a_r));


%% Initializing Lattice Conditions


i_alpha=80;%150,110,70
if i_alpha<42
    pol = ['P=0'];
elseif i_alpha<62
    pol = ['P=a_2'];
elseif i_alpha<121
    pol = ['P=a_2-a_1'];
else
    pol = ['P=0'];
end

%choose the same angles from the homogeneous family when origin is at the hard edge.
alpha0=Alpha(i_alpha);
gamma0=Gamma(i_alpha);
theta0=Theta(i_alpha);

%Convert them into the initial angles solved fot the origin at the soft edge.
alpha_0=2*pi-alpha0-psi_cb-psi_ar;
theta_0=2*pi-theta0-psi_bb-psi_cr;
gamma_0=2*pi-gamma0-psi_ab-psi_br;

%Create angle matrice
theta_matrix=zeros(n+1,m);
alpha_matrix=zeros(n+1,m);
gamma_matrix=zeros(n+1,m);
kappa_matrix=zeros(n+1,m);
phi_matrix=zeros(n+1,m);
%Store angles of the complementary angles
theta_matrix0=zeros(n+1,m);
alpha_matrix0=zeros(n+1,m);
gamma_matrix0=zeros(n+1,m);

theta_matrix(2,:)=theta_0;
theta_matrix(2:n+1,1)=theta_0;
% theta_matrix(2,1)=theta_0;
alpha_matrix(3:n+1,1)=alpha_0+0.1;
% alpha_matrix(3,1)=alpha_0;
gamma_matrix(2,2:m)=gamma_0;

x = linspace(0,2*pi,m);
k = 2; %spatial wavenumber
wantplot = 0;
cut_ind = m; %index to cut on right side
cut_indy =0; %index to cut y on top (m-cut_indy)

ep =0* 1e-3;
wantvid = 0; %one if want video 0 o.w.
% pert= ep*rand(1,m);
pert = ep*sin(k.*x);
ky_res = 2*pi/(n);
ky_title = [' ky= [-0.041-0.0004i, -0.251-0.046i]'];
ky_realfit = -[-0.0411216 -0.250518]; %kx=0.0524
ky_imfit = -[-0.000412539 -0.0460207];
%ky_title = [' ky= [-0.08-0.001i, -0.452-0.16i]'];
%ky_realfit = [-0.0822992 -0.451516];
%ky_imfit = [-0.0016498 -0.160829];
kx_in = 2*pi/(m/k);
% Give initial conditions at the left and bottom boundaries
pol_mat = zeros(n,m);

theta_matrix(2,:)=theta_matrix(2,:)+pert;
%theta_matrix(2,m) = theta_0+1e-6;

%% Solve hexagon coordinates and angles of each unit cell
for i=2:n
    while 1
    for j=1:m-1
        [hex_coord,theta_matrix(i+1,j+1),alpha_matrix(i+1,j+1),...
            gamma_matrix(i+1,j+1),phi_matrix(i+1,j+1)]=...
            solve_hexagon_analytical(theta_matrix(i,j),alpha_matrix(i+1,j),gamma_matrix(i,j+1));
        %         [hex_coord,theta_matrix(i+1,j+1),alpha_matrix(i+1,j+1),...
        %             gamma_matrix(i+1,j+1),phi_matrix(i+1,j+1)]=...
        %             solve_hexagon_Newton(theta_matrix(i,j),alpha_matrix(i+1,j),gamma_matrix(i,j+1));
        Hexagon_coord(i,j,:)=hex_coord;
    end
    if abs(alpha_matrix(i+1,1)-alpha_matrix(i+1,m))<1e-8
        break
    end
    %for fixed conditions comment these two out and set the left boundary
%     theta_matrix(i+1,1)=theta_matrix(i,j+1);
    alpha_matrix(i+1,1)=alpha_matrix(i+1,m); 
    end
end





%Solve kappa and phi of the bottom edge
i=2;
for j=1:m-1
    phi_matrix(i,j+1)=gamma_matrix(i,j+1)+theta_matrix(i,j+1)+psi_ab+psi_bb-2*pi;
    kappa_matrix(i,j+1)=psi_ab+gamma_matrix(i,j+1)-pi;
end


%Solve kappa and phi of the left edge
j=1;
for i=2:n
    phi_matrix(i+1,j)=2*pi-alpha_matrix(i+1,j)-theta_matrix(i,j)-psi_ar-psi_cr;
    kappa_matrix(i+1,j)=3*pi-alpha_matrix(i+1,j)-theta_matrix(i,j)-psi_ar-psi_cr-psi_bb;
end





theta_matrix0(2:n+1,:)=2*pi-theta_matrix(2:n+1,:)-psi_bb-psi_cr;
theta_matrix0=rot90(theta_matrix0,2);
% theta_matrix0(n+1,:)=theta_matrix0(n,:);

alpha_matrix0(3:n+1,:)=2*pi-alpha_matrix(3:n+1,:)-psi_cb-psi_ar;
alpha_matrix0=rot90(alpha_matrix0,2);
% alpha_matrix0(n+1,:)=alpha_matrix0(n,:);

gamma_matrix0(2:n+1,2:m)=2*pi-gamma_matrix(2:n+1,2:m)-psi_ab-psi_br;
gamma_matrix0=rot90(gamma_matrix0,2);
% gamma_matrix0(n+1,:)=gamma_matrix0(n,:);


for i = 1:n
    for j = 1:m
        if theta_matrix0(i,j)<=Theta(122)
            pol_mat(i,j) = 0;
        elseif theta_matrix0(i,j)<=Theta(63)
            pol_mat(i,j) = 1;
        elseif theta_matrix0(i,j) <= Theta(44)
            pol_mat(i,j) = 0.5;
        else
            pol_mat(i,j) = 0;
        end
        
    end
end

save('a2a1_fixed_bc.mat')