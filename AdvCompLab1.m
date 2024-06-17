%Q1 axial strain in x direction
clear;
load US_data.dat;
axial_strain_x = zeros(51,46); %make 0 matrix for axial strain
h = 0.001; %x direction interval for calculation

for x = 2:47 % make vector with indexes of where x changes
    z(1,1) = 0;
    z(x,1) = 51*(x-1);
end

for y = 1:46 % split up sections of data for easier calculation
    organizedMatrix{y} = US_data(z(y,1)+1:z(y+1,1),:);
end

for j = 1:51 %rows
    for i = 1:46 %columns
        if i == 1
            axial_strain_x(j,i) = (organizedMatrix{i+1}(j,4) - organizedMatrix{i}(j,4))/h;
        elseif i == 46
            axial_strain_x(j,i) = (organizedMatrix{i}(j,4) - organizedMatrix{i-1}(j,4))/h;
        else
            axial_strain_x(j,i) = (organizedMatrix{i+1}(j,4) - organizedMatrix{i-1}(j,4))/(2*h); %math on each cell not on BC
        end
    end
end

minrow = min(min(axial_strain_x(:,1:46)));
maxrow = max(max(axial_strain_x(:,1:46)));

figure;imagesc(axial_strain_x);caxis([minrow maxrow]);title('x axial strain');colorbar

%% axial strain in y direction
clear;
load US_data.dat;
axial_strain_y = zeros(51,46); %make 0 matrix for y axial strain
h = 0.001;

for p = 1:51 % group y values together
    counter = 1;
    for x = p:51:2346 % populate organized matrix in y direction
    organizedMatrix_y{p,1}(counter,:) = US_data(x,:);
    counter = counter + 1;
    end
end

for j = 1:51 % all columns
    for i = 1:46 % rows
        if j == 1
            axial_strain_y(j,i) = (organizedMatrix_y{j+1,1}(i,5) - organizedMatrix_y{j,1}(i,5))/h;
        elseif j == 51
            axial_strain_y(j,i) = (organizedMatrix_y{j,1}(i,5) - organizedMatrix_y{j-1,1}(i,5))/h;
        else
            axial_strain_y(j,i) = (organizedMatrix_y{j+1,1}(i,5) - organizedMatrix_y{j-1,1}(i,5))/(2*h);
        end
    end
end

mincol = min(min(axial_strain_y(:,1:46)));
maxcol = max(max(axial_strain_y(:,1:46)));

figure;imagesc(axial_strain_y);caxis([mincol maxcol]);title('y axial strain');colorbar

%% shear strain
clear;
load US_data.dat;
axial_strain_x = zeros(51,46); %make 0 matrix for axial strain
h = 0.001; %x direction interval for calculation

for x = 2:47 % make vector with indexes of where x changes
    z(1,1) = 0;
    z(x,1) = 51*(x-1);
end

for y = 1:46 % split up sections of data for easier calculation
    organizedMatrix{y} = US_data(z(y,1)+1:z(y+1,1),:);
end

for j = 1:51 %rows
    for i = 1:46
        if i == 1
            axial_strain_x(j,i) = (organizedMatrix{i+1}(j,5) - organizedMatrix{i}(j,5))/h;
        elseif i == 46
            axial_strain_x(j,i) = (organizedMatrix{i}(j,5) - organizedMatrix{i-1}(j,5))/h;
        else
            axial_strain_x(j,i) = (organizedMatrix{i+1}(j,5) - organizedMatrix{i-1}(j,5))/(2*h); %math on each cell not on BC
        end
    end
end

axial_strain_y = zeros(51,46); %make 0 matrix for y axial strain
h = 0.001;

for p = 1:51 % group y values together
    counter = 1;
    for x = p:51:2346 % populate organized matrix in y direction
    organizedMatrix_y{p,1}(counter,:) = US_data(x,:);
    counter = counter + 1;
    end
end

for j = 1:51 % all columns
    for i = 1:46 % rows
        if j == 1
            axial_strain_y(j,i) = (organizedMatrix_y{j+1,1}(i,4) - organizedMatrix_y{j,1}(i,4))/h;
        elseif j == 51
            axial_strain_y(j,i) = (organizedMatrix_y{j,1}(i,4) - organizedMatrix_y{j-1,1}(i,4))/h;
        else
            axial_strain_y(j,i) = (organizedMatrix_y{j+1,1}(i,4) - organizedMatrix_y{j-1,1}(i,4))/(2*h);
        end
    end
end

shear_strain = (axial_strain_y + axial_strain_x)/2;

minshear = min(min(shear_strain));
maxshear = max(max(shear_strain));
figure;imagesc(shear_strain);caxis([minshear maxshear]);title('shear strain');colorbar

%% Q3 Jacobi
clear;
U=zeros(21,21); %make 0 matrix
error = 1;
itr = 0;
for x = 1:21 % right wall BC
    U(x,21) = 1;
end
for y = 1:21 % top and bottom BC
    U(1,y) = 1;
    U(21,y) = 1;
end
for z = 0.05:0.05:0.95 % left wall BC
    z_index = int64(z * 20 + 1);
    U(z_index,1) = cos(2*pi*z);
end

while (error > 1e-5 && itr < 10000)
    itr=itr+1;
    Uold=U;
    for i=0.05:0.05:.95 %
        index_i = int64(i*20 + 1); % for indexing
        for j=0.05:0.05:.95
            index_j = int64(j*20 + 1); % for indexing
            U(index_i,index_j)=1/4*(Uold(index_i-1,index_j)+Uold(index_i+1,index_j)+Uold(index_i,index_j-1)+Uold(index_i,index_j+1));
        end
    end
    errorold=error;
    error=max(max(abs(U-Uold)));
end
fprintf('Point Iterative - Jacobi %d\n',itr);
fprintf('Value at 0.7,0.7 - %d\n',U(15,15));
figure;imagesc(U);caxis([0 1]);title('Jacobi');colorbar


%% Gauss-Seidel
clear;
U=zeros(21,21); %make 0 matrix
for x = 1:21 % right wall BC
    U(x,21) = 1;
end
for y = 1:21 % top and bottom BC
    U(1,y) = 1;
    U(21,y) = 1;
end
for z = 0.05:0.05:0.95 % left wall BC
    z_index = int64(z * 20 + 1);
    U(z_index,1) = cos(2*pi*z);
end

error = 1;
itr = 0;
while (error > 1e-5 && itr < 10000)
    itr=itr+1;
    Uold=U;
    for i=0.05:0.05:.95
        index_i = int64(i*20 + 1); % for indexing
        for j=0.05:0.05:.95
            index_j = int64(j*20 + 1); % for indexing
            U(index_i,index_j)=1/4*(U(index_i-1,index_j)+Uold(index_i+1,index_j)+U(index_i,index_j-1)+Uold(index_i,index_j+1));
        end
    end
    errorold=error;
    error=max(max(abs(U-Uold)));
end
fprintf('Point Iterative - Gauss-Seidel %d\n',itr);
fprintf('Value at 0.7,0.7 - %d\n',U(15,15));
figure;imagesc(U);caxis([0 1]);title('Gauss-Seidel');colorbar

%% SOR method
clear;
h = .05;
a = 1.05;
b = 1.05;

rowGauss = 1/2*(cos(pi*h/a) + cos(pi*h/b));
w = 2/(1 + sqrt(1 - (rowGauss)^2));


U=zeros(21,21); %make 0 matrix
for x = 1:21 % right wall BC
    U(x,21) = 1;
end
for y = 1:21 % top and bottom BC
    U(1,y) = 1;
    U(21,y) = 1;
end
for z = 0.05:0.05:0.95 % left wall BC
    z_index = int64(z * 20 + 1);
    U(z_index,1) = cos(2*pi*z);
end

error = 1;
itr = 0;
while (error > 1e-5 && itr < 10000)
    itr=itr+1;
    Uold=U;
    for i=0.05:0.05:.95
        index_i = int64(i*20 + 1); % for indexing
        for j=0.05:0.05:.95
            index_j = int64(j*20 + 1); % for indexing
            U(index_i,index_j)=w/4*(U(index_i-1,index_j)+Uold(index_i+1,index_j)+U(index_i,index_j-1)+Uold(index_i,index_j+1))+((1-w)*Uold(index_i,index_j));
        end
    end
    errorold=error;
    error=max(max(abs(U-Uold)));
end
fprintf('Point Iterative - SOR %d\n',itr);
fprintf('Value at 0.7,0.7 - %d\n',U(15,15));
figure;imagesc(U);caxis([0 1]);title('SOR');colorbar

wPlot=zeros(2001,2);
for z = 0:.001:2
wPlot(int64(z*1000+1),1) = z;
end

for w = 0:0.001:2
    U=zeros(21,21); %make 0 matrix
for x = 1:21 % right wall BC
    U(x,21) = 1;
end
for y = 1:21 % top and bottom BC
    U(1,y) = 1;
    U(21,y) = 1;
end
for z = 0.05:0.05:0.95 % left wall BC
    z_index = int64(z * 20 + 1);
    U(z_index,1) = cos(2*pi*z);
end
error = 1;
itr = 0;
while (error > 1e-5 && itr < 10000)
    itr=itr+1;
    Uold=U;
    for i=0.05:0.05:.95
        index_i = int64(i*20 + 1); % for indexing
        for j=0.05:0.05:.95
            index_j = int64(j*20 + 1); % for indexing
            U(index_i,index_j)=w/4*(U(index_i-1,index_j)+Uold(index_i+1,index_j)+U(index_i,index_j-1)+Uold(index_i,index_j+1))+((1-w)*Uold(index_i,index_j));
        end
    end
    errorold=error;
    error=max(max(abs(U-Uold)));
end
wPlot(int64(w*1000 + 1),2) = itr;
end

plot(wPlot(:,1),wPlot(:,2));
xlabel('w value')
ylabel('Iterations')
grid on
title('Total iterations at w values from 0-2')
