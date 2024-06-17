counter1 = 1;
error = 1;
t = 0.004;
h = 0.006;
D = 0.001;
itr = 0;
i = 1;

for x = 0:.006:.3 % t=0 intital condition
    U(1,counter1) = (2.5/(sqrt(2*pi))*exp(-((30*x-4.5)^2)/2));
    counter1 = counter1 + 1;
end

while (error > 1e-4 && itr < 10000)
    itr = itr + 1;   
    i = i + 1;
    for j = 1:51
        if j == 1
            U(i,j) = D*t*((2*U(i-1,j+1)-2*U(i-1,j))/h^2)+U(i-1,j);
        elseif j == 51
            U(i,j) = D*t*((-2*U(i-1,j)+2*U(i-1,j-1))/h^2)+U(i-1,j);
        else
            U(i,j) = D*t*((U(i-1,j+1)-2*U(i-1,j)+U(i-1,j-1))/h^2)+U(i-1,j);
        end
    end
    errorold=error;
    error=(max(abs(U(i,:)-U(i-1,:))))/(max(abs(U(i,:))));
end

axis = 0:.006:.3;

plot(axis,U(1,:))
hold on
plot(axis,U(51,:))
hold on
plot(axis,U(101,:))
hold on
plot(axis,U(151,:))
hold on
plot(axis,U(251,:))
hold on
plot(axis,U(451,:))
hold on
plot(axis,U(751,:))
hold on
plot(axis,U(1001,:))
hold on
plot(axis,U(1857,:))
ylabel('Concentration')
xlabel('x position')
legend('t=0','t=0.2','t=0.4','t=0.6','t=1.0','t=1.8','t=3.0','t=4.0','t=7.436')
grid on

%check threshold

plot(axis,U(100,:)) %fxn breaks after t > 0.018

%% 
load LINE_LIST_1.DAT
load POINT_LIST_1.DAT

line = LINE_LIST_1;
point = POINT_LIST_1;
total = 0;
length = 0;
answer = 0;

for i = 1:464
    node_ind = find(line(:,2) == i); % find index in line for node of interest (NOI)
    node_minus1 = line(node_ind,3); % find index in point for NOI - 1
    plus_node_ind = find(line(:,3) == i); % find index in line for NOI + 1
    node_plus1 = line(plus_node_ind,2); % find index in point for NOI + 1
    n_vectorx = point(node_minus1,2) - point(node_plus1,2); % final - initial for x position
    n_vectory = point(node_minus1,3) - point(node_plus1,3); % final - initial for y position
    n_result = [n_vectorx n_vectory 0];
    z_vector = [0 0 1];
    n1 = cross(n_result,z_vector);
    n1 = n1/norm(n1); % normalized vector
    j1 = [point(i,4) point(i,5) 0]; % flux vector
    total1 = dot(j1,n1);

    % find second half of equation
    node_doubleCCW_ind = find(line(:,3 == node_minus1));
    node_doubleCCW = line(node_doubleCCW_ind,2);
    
    n2_vectorx = point(node_doubleCCW,2) - point(i,2); % final - initial for x position
    n2_vectory = point(node_doubleCCW,3) - point(i,3); % final - initial for y position
    n2_result = [n2_vectorx n2_vectory 0];
    z2_vector = [0 0 1];
    n2 = cross(n2_result,z2_vector);
    n2 = n2/norm(n2); % normalized vector
    j2 = [point(node_minus1,4) point(node_minus1,5) 0]; % flux vector
    total2 = dot(j2,n2);

    total = total1 + total2;

    % calculate the length of the entire figure
    current_node = [point(i,2) point(i,3)];
    current_node_m1 = [point(node_minus1,2) point(node_minus1,3)];
    length = norm(current_node_m1-current_node);
    answer = answer + total*(length/2);
end 



%%
load LINE_LIST_2.DAT
load POINT_LIST_2.DAT

line = LINE_LIST_2;
point = POINT_LIST_2;
total = 0;
length = 0;

for i = 1:464
    node_ind = find(line(:,2) == i); % find index in line for node of interest
    node_minus1 = line(node_ind,3); % find index for point
    plus_node_ind = find(line(:,3) == i);
    node_plus1 = line(plus_node_ind,2);
    n_vectorx = point(node_minus1,2) - point(node_plus1,2);
    n_vectory = point(node_minus1,3) - point(node_plus1,3);
    n_result = [n_vectorx n_vectory 0];
    z_vector = [0 0 1];
    n1 = cross(n_result,z_vector);
    j1 = [point(i,4) point(i,5) 0];
    total = total + dot(j1,n1);

    % calculate the length of the entire figure
    current_node = [point(i,2) point(i,3)];
    current_node_m1 = [point(node_minus1,2) point(node_minus1,3)];
    length = length + norm(current_node_m1-current_node);
end 

answer = total*(length/2);






