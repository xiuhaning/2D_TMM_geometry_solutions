clear all
clc
load a2a1_soft_20rows.mat

%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons

% m1 = n;
% n1 = m;
% m = m1;
% n = n1;
% m = 10;
% n = 10;
Coor_hexagon_x=zeros(n,m-1,7); 
Coor_hexagon_y=zeros(n,m-1,7);
% x, y coordinates for each hexagon from A, B, C, D, E, F, to A 
Coor_unit_cell_x=zeros(n+1,m,9);
Coor_unit_cell_y=zeros(n+1,m,9);
% x, y coordinates for each unit cell from C, A, B, C, A', B' C, P, to Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blue Triangle
global a_b; global b_b; global c_b;
% a_b=0.5;b_b=0.7;c_b=1;
a_b=0.8;b_b=1;c_b=0.4;

global psi_ab; global psi_bb; global psi_cb;
psi_ab=acos((a_b^2-b_b^2-c_b^2)/(-2*b_b*c_b));
psi_bb=acos((b_b^2-a_b^2-c_b^2)/(-2*a_b*c_b));
psi_cb=acos((c_b^2-b_b^2-a_b^2)/(-2*b_b*a_b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Red Triangle
global a_r; global b_r; global c_r;
% a_r=0.4;b_r=0.8;c_r=1;
a_r=1;b_r=0.5;c_r=0.7;

global psi_ar; global psi_br; global psi_cr;
psi_ar=acos((a_r^2-b_r^2-c_r^2)/(-2*b_r*c_r));
psi_br=acos((b_r^2-a_r^2-c_r^2)/(-2*a_r*c_r));
psi_cr=acos((c_r^2-b_r^2-a_r^2)/(-2*b_r*a_r));

% theta_c1=pi-psi_ar
% theta_c3=pi-psi_cb
% theta_c2=psi_br+psi_bb

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Represent an unit cell at the origin: theta(2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=2;j=1;
A = [0,0];
B = [a_r*cos(psi_bb+theta_matrix(i,j)),-a_r*sin(psi_bb+theta_matrix(i,j))];
C = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr), -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)]; 
Aa = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*cos(3*pi-alpha_matrix(i+1,j)...
    -theta_matrix(i,j)-psi_ar-psi_cr-psi_bb),...
    -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*sin(3*pi-alpha_matrix(i+1,j)-...
    theta_matrix(i,j)-psi_ar-psi_cr-psi_bb)];
Bb = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)-...
    b_b*cos(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar),...
    -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+...
    b_b*sin(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar)];
P=(A+C)/2;
Q=(Aa+C)/2;
Coor_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]';
Coor_unit_cell_x(i,j,:)=Coor_unit_cell(1,:);
Coor_unit_cell_y(i,j,:)=Coor_unit_cell(2,:);


C = [a_b*cos(psi_bb), -a_b*sin(psi_bb)]; 
A = [0,0];
B = [c_b,0];
Coor_unit_cell=[0,0;0,0;0,0;C;A;B;C;0,0;0,0]';
Coor_unit_cell_x(i-1,j,:)=Coor_unit_cell(1,:);
Coor_unit_cell_y(i-1,j,:)=Coor_unit_cell(2,:);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce unit cells of the 1st and 2nd row of lattice (bottom edge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rotation=[1,0;0,1];
hexcoor=0;
i=2;
for j=2:m 
    if theta_matrix(i,j)==0
        break;
    end
    %2nd row
    A = [0,0];
    B = [a_r*cos(psi_bb+theta_matrix(i,j)),-a_r*sin(psi_bb+theta_matrix(i,j))];
    C = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr), -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)]; 
    Aa = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*cos(3*pi-alpha_matrix(i+1,j)...
        -theta_matrix(i,j)-psi_ar-psi_cr-psi_bb),...
        -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*sin(3*pi-alpha_matrix(i+1,j)-...
        theta_matrix(i,j)-psi_ar-psi_cr-psi_bb)];
    Bb = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)-...
    b_b*cos(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar),...
    -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+...
    b_b*sin(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar)];
    P=(A+C)/2;
    Q=(Aa+C)/2;
    Coor_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]';
    
    %rotate from theta_i,j to theta_i,j+1
    %rotation angle should be repeatedly added ever loop
    rotation1=[cos(phi_matrix(i,j)),-sin(phi_matrix(i,j));
         sin(phi_matrix(i,j)),cos(phi_matrix(i,j))]*rotation;
    %the coordinates of theta(i,j) is rotated back to the origion system (thet(2,1))
    Coor_unit_cell=rotation1*Coor_unit_cell; 
    % find the distance of the origin of the local system to the origin of
    % theta(2,1)
    hexcoor=hexcoor+rotation*...
        [Hexagon_coord(i,j-1,3)+a_r*cos(kappa_matrix(i,j));Hexagon_coord(i,j-1,4)+a_r*sin(kappa_matrix(i,j))];
    Coor_unit_cell=Coor_unit_cell+hexcoor;

    Coor_unit_cell_x(i,j,:)=Coor_unit_cell(1,:);
    Coor_unit_cell_y(i,j,:)=Coor_unit_cell(2,:);
    
    %1st row, only blue triangles
    C = [a_b*cos(psi_bb), -a_b*sin(psi_bb)]; 
    A = [0,0];
    B = [c_b,0];
    Coor_unit_cell=[0,0;0,0;0,0;C;A;B;C;0,0;0,0]';
    Coor_unit_cell=rotation1*Coor_unit_cell; 
    Coor_unit_cell=Coor_unit_cell+hexcoor;
    
    Coor_unit_cell_x(i-1,j,:)=Coor_unit_cell(1,:);
    Coor_unit_cell_y(i-1,j,:)=Coor_unit_cell(2,:);
    
    rotation=rotation1;   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce unit cells of the 1st column of lattice (left edge)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rotation=[1,0;0,1];
hexcoor=0;
j=1;
for i=3:n+1
    if theta_matrix(i,j)==0
        break;
    end
    A = [0,0];
    B = [a_r*cos(psi_bb+theta_matrix(i,j)),-a_r*sin(psi_bb+theta_matrix(i,j))];
    C = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr), -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)]; 
    if i<n+1 %the n+1 row only has a red triangle
    Aa = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*cos(3*pi-alpha_matrix(i+1,j)...
        -theta_matrix(i,j)-psi_ar-psi_cr-psi_bb),...
        -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*sin(3*pi-alpha_matrix(i+1,j)-...
        theta_matrix(i,j)-psi_ar-psi_cr-psi_bb)];
    Bb = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)-...
        b_b*cos(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar),...
        -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+...
        b_b*sin(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar)];
    P=(A+C)/2;
    Q=(Aa+C)/2;
    end
    Coor_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]';
    
    rotation1=[cos(phi_matrix(i,j)),-sin(phi_matrix(i,j));
        sin(phi_matrix(i,j)),cos(phi_matrix(i,j))]*rotation;
    %the coordinates of theta(i,j) is rotated back to the origion system (thet(2,1))
    Coor_unit_cell=rotation1*Coor_unit_cell; 
    % find the distance of the origin of the local system to the origin of
    % theta(2,1)
    hexcoor=hexcoor+rotation*...
        [Hexagon_coord(i-1,j,11)+a_b*cos(kappa_matrix(i,j));Hexagon_coord(i-1,j,12)+a_b*sin(kappa_matrix(i,j))];
    Coor_unit_cell=Coor_unit_cell+hexcoor;
    if i<n+1
        Coor_unit_cell_x(i,j,:)=Coor_unit_cell(1,:);
        Coor_unit_cell_y(i,j,:)=Coor_unit_cell(2,:);
    else %n+1 row, only red triangles
        Coor_unit_cell_x(i,j,1:4)=Coor_unit_cell(1,1:4);
        Coor_unit_cell_y(i,j,1:4)=Coor_unit_cell(2,1:4);
    end
    
    rotation=rotation1;   
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce unit cells diagonally from the bottom edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotation_bottom=[1,0;0,1];
for k=2:m %generate unit cells diagonally from the bottom edge
    i=3;j=k;
    rotation_bottom=[cos(phi_matrix(2,k-1)),-sin(phi_matrix(2,k-1));
            sin(phi_matrix(2,k-1)),cos(phi_matrix(2,k-1))]*rotation_bottom;
    rotation=rotation_bottom; %rotate the angles from the bottom edge (theta(2,k)) to the origin
    hexcoor=0;
    while i<=n+1 && j<=m
        if theta_matrix(i,j)==0
            break;  
        end
        A = [0,0];
        B = [a_r*cos(psi_bb+theta_matrix(i,j)),-a_r*sin(psi_bb+theta_matrix(i,j))];
        C = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr), -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)]; 
        if i<n+1 %the n+1 row only has a red triangle
            Aa = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*cos(3*pi-alpha_matrix(i+1,j)...
                -theta_matrix(i,j)-psi_ar-psi_cr-psi_bb),...
                -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*sin(3*pi-alpha_matrix(i+1,j)-...
                theta_matrix(i,j)-psi_ar-psi_cr-psi_bb)];
            Bb = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)-...
                b_b*cos(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar),...
                -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+...
                b_b*sin(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar)];
            P=(A+C)/2;
            Q=(Aa+C)/2;
        end
        Coor_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]';
    
        rotation1=[cos(phi_matrix(i,j)),-sin(phi_matrix(i,j));
            sin(phi_matrix(i,j)),cos(phi_matrix(i,j))]*rotation;
        %the coordinates of theta(i,j) is rotated back to the origion system (thet(2,1))
        Coor_unit_cell=rotation1*Coor_unit_cell; 
        % find the distance of the origin of the local system to the origin of
        % theta(2,1)
        hexcoor=hexcoor+rotation*[Hexagon_coord(i-1,j-1,7);Hexagon_coord(i-1,j-1,8)];
        Coor_unit_cell=Coor_unit_cell+hexcoor;

        if i<n+1
            Coor_unit_cell_x(i,j,:)=Coor_unit_cell(1,:);
            Coor_unit_cell_y(i,j,:)=Coor_unit_cell(2,:);
        else
            Coor_unit_cell_x(i,j,1:4)=Coor_unit_cell(1,1:4);
            Coor_unit_cell_y(i,j,1:4)=Coor_unit_cell(2,1:4);
        end
        Coor_unit_cell_x(i,j,:)=Coor_unit_cell_x(i,j,:)+Coor_unit_cell_x(2,k-1,2);
        Coor_unit_cell_y(i,j,:)=Coor_unit_cell_y(i,j,:)+Coor_unit_cell_y(2,k-1,2);
        i=i+1;j=j+1;
        rotation=rotation1; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce unit cells diagonally from the left edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotation_left=[1,0;0,1];
for k=3:n+1 %generate unit cells diagonally from the bottom edge
    i=k;j=2;
    rotation_left=[cos(phi_matrix(k-1,1)),-sin(phi_matrix(k-1,1));
            sin(phi_matrix(k-1,1)),cos(phi_matrix(k-1,1))]*rotation_left;
    rotation=rotation_left; %rotate the angles from the left edge (theta(k,1)) to the origin
    hexcoor=0;
    while i<=n+1 && j<=m
        if theta_matrix(i,j)==0
            break;  
        end
        A = [0,0];
        B = [a_r*cos(psi_bb+theta_matrix(i,j)),-a_r*sin(psi_bb+theta_matrix(i,j))];
        C = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr), -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)]; 
        if i<n+1 %the n+1 row only has a red triangle
            Aa = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*cos(3*pi-alpha_matrix(i+1,j)...
                -theta_matrix(i,j)-psi_ar-psi_cr-psi_bb),...
                -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+a_b*sin(3*pi-alpha_matrix(i+1,j)-...
                theta_matrix(i,j)-psi_ar-psi_cr-psi_bb)];
            Bb = [b_r*cos(theta_matrix(i,j)+psi_bb+psi_cr)-...
                b_b*cos(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar),...
                -b_r*sin(theta_matrix(i,j)+psi_bb+psi_cr)+...
                b_b*sin(theta_matrix(i,j)+psi_bb+psi_cr+alpha_matrix(i+1,j)+psi_cb+psi_ar)];
            P=(A+C)/2;
            Q=(Aa+C)/2;
        end
        Coor_unit_cell=[C;A;B;C;Aa;Bb;C;P;Q]';
    
        rotation1=[cos(phi_matrix(i,j)),-sin(phi_matrix(i,j));
            sin(phi_matrix(i,j)),cos(phi_matrix(i,j))]*rotation;
        %the coordinates of theta(i,j) is rotated back to the origion system (thet(2,1))
        Coor_unit_cell=rotation1*Coor_unit_cell; 
        % find the distance of the origin of the local system to the origin of
        % theta(2,1)
        hexcoor=hexcoor+rotation*[Hexagon_coord(i-1,j-1,7);Hexagon_coord(i-1,j-1,8)];
        Coor_unit_cell=Coor_unit_cell+hexcoor;

        if i<n+1
            Coor_unit_cell_x(i,j,:)=Coor_unit_cell(1,:);
            Coor_unit_cell_y(i,j,:)=Coor_unit_cell(2,:);
        else
            Coor_unit_cell_x(i,j,1:4)=Coor_unit_cell(1,1:4);
            Coor_unit_cell_y(i,j,1:4)=Coor_unit_cell(2,1:4);
        end
        Coor_unit_cell_x(i,j,:)=Coor_unit_cell_x(i,j,:)+Coor_unit_cell_x(k-1,1,2);
        Coor_unit_cell_y(i,j,:)=Coor_unit_cell_y(i,j,:)+Coor_unit_cell_y(k-1,1,2);
        i=i+1;j=j+1;
        rotation=rotation1; 
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reproduce hexagons from left bottom to right top (2,1) to (n,m-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:n
    for j=1:m-1
        Coor_hexagon_x(i,j,:)=[Coor_unit_cell_x(i,j,2),Coor_unit_cell_x(i,j+1,3),...
            Coor_unit_cell_x(i,j+1,1),Coor_unit_cell_x(i+1,j+1,2),Coor_unit_cell_x(i+1,j+1,3),...
            Coor_unit_cell_x(i,j,1),Coor_unit_cell_x(i,j,2)]; 
        Coor_hexagon_y(i,j,:)=[Coor_unit_cell_y(i,j,2),Coor_unit_cell_y(i,j+1,3),...
            Coor_unit_cell_y(i,j+1,1),Coor_unit_cell_y(i+1,j+1,2),Coor_unit_cell_y(i+1,j+1,3),...
            Coor_unit_cell_y(i,j,1),Coor_unit_cell_y(i,j,2)]; 
    end
end


%%
%Plot all unit cells and hexagons to form a n-by-m lattice
rotation180=[cos(pi),-sin(pi);
            sin(pi),cos(pi)];

figure;
for i=1:n+1
    for j=1:m
        for k=1:9
            XY_unit(1,k)=Coor_unit_cell_x(i,j,k);XY_unit(2,k)=Coor_unit_cell_y(i,j,k);
        end
%         XY_unit=rotation180*XY_unit;
        if i==1
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
            hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.0);
        elseif i==n+1           
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
            hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.0);
        else
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
            hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.0);
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
            hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.0);
%             hold on;plot(XY_unit(1,8:9),XY_unit(2,8:9),'g-','linewidth',2)
        end
    end
end
axis equal
        
        
        
for i=1:n
    for j=1:m
        for k=1:9
            XY_unit(1,k)=Coor_unit_cell_x(i,j,k);XY_unit(2,k)=Coor_unit_cell_y(i,j,k);
        end
        XY_unit=rotation180*XY_unit;
        for k=1:9
            Coor_unit_cell_x_new(i,j,k)=XY_unit(1,k);Coor_unit_cell_y_new(i,j,k)=XY_unit(2,k);
        end
%         if i==1
% %             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
%             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.5);
%         elseif i==n+1           
% %             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
%             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.5);
%         else
% %             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
%             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.5);
% %             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
%             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.5);
% %             hold on;plot(XY_unit(1,8:9),XY_unit(2,8:9),'g-','linewidth',2)
%         end
    end
end

 
%  figure;
% for i=1:n+1
%     for j=1:m
%         for k=1:9
%             XY_unit(1,k)=Coor_unit_cell_x_new(i,j,k);XY_unit(2,k)=Coor_unit_cell_y_new(i,j,k);
%         end
%         if i==1
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
% %             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.5);
%         elseif i==n+1           
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
% %             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.5);
%         else
%             hold on;fill(XY_unit(1,1:4),XY_unit(2,1:4),'b')
% %             hold on;plot(XY_unit(1,1:4),XY_unit(2,1:4),'b-','linewidth',1.5);
%             hold on;fill(XY_unit(1,4:7),XY_unit(2,4:7),'r')
% %             hold on;plot(XY_unit(1,4:7),XY_unit(2,4:7),'r-','linewidth',1.5);
% %             hold on;plot(XY_unit(1,8:9),XY_unit(2,8:9),'g-','linewidth',2)
%         end
%     end
% end
%  axis equal

%%

% figure;
% pert_mat1 = (alpha_matrix0(17:30,1:cut_ind)-alpha0);
% pert_mat1 = rot90(pert_mat1,2);
% for i=2:13
%     for j=1:m-1
%         for k=1:7
%             XY_hex(1,k)=Coor_hexagon_x(i,j,k);XY_hex(2,k)=Coor_hexagon_y(i,j,k);
%         end
%         %scatter(-Coor_hexagon_x(i,j,1),-Coor_hexagon_y(i,j,1),50,pert_mat(i,j),'filled')
%         XY_hex=rotation180*XY_hex;
%         hold on; plot(XY_hex(1,:),XY_hex(2,:),'k-','linewidth',0.4,'MarkerSize',0.5)
%         x_coord = squeeze(Coor_hexagon_x(i,j,:));
%         y_coord = squeeze(Coor_hexagon_y(i,j,:));
% %         hold on;fill(x_coord,y_coord,pert_mat1(i,j))
%         hold on;fill(XY_hex(1,:),XY_hex(2,:),pert_mat1(i,j))
%     end
% end
% axis equal
% x_q = Coor_hexagon_x(n,floor(m/2));
% y_q = Coor_hexagon_y(n,floor(m/2));
% u_q1 = Coor_hexagon_x(n-1,floor(m/2)-1)-x_q;
% u_q2 = Coor_hexagon_x(n-1,floor(m/2)-1)-x_q;
% u_q3 = Coor_hexagon_x(n,floor(m/2)+1)-x_q;
% U = -[u_q1,u_q1+u_q3,-u_q3];
% v_q1 = Coor_hexagon_y(n-1,floor(m/2))-y_q;
% v_q2 = Coor_hexagon_y(n-1,floor(m/2)+1)-y_q;
% v_q3 = Coor_hexagon_y(n,floor(m/2)+1)-y_q;
% V = -[v_q1, v_q1+v_q3, -v_q3];
% quiver(-[x_q x_q x_q], -[y_q y_q y_q], (m/5)*U, (n/5)*V,'LineWidth',2,'Color','m') 
%title('Lattice with abs(\theta pert)')
% colorbar
% colormap jet
% shading flat
% axis square



