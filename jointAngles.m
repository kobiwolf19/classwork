function [a,b,y] = jointAngles(prox,dist)

angles = prox'*dist;

a = atan2(angles(2,1), angles(1,1));
a = rad2deg(a);
b = atan2(-angles(3,1), sqrt((angles(1,1)^2 + angles(2,1)^2)));
b = rad2deg(b);
y = atan2(angles(3,2), angles(3,3));
y = rad2deg(y);
end
