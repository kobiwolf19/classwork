function [K] = HelicalAngles(prox,dist)

angles = prox'*dist;

theta = acos((angles(1,1) + angles(2,2) + angles(3,3) - 1)/2);

r = [angles(3,2) - angles(2,3); angles(1,3) - angles(3,1); angles(2,1) - angles(1,2)];

K = (1/(2*sin(theta))) * r;
theta = rad2deg(theta);

K = theta*K;

end
