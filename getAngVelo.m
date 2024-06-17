function [angVelo,angAccel] = getAngVelo(R)

for i = 1:length(R)
    L1(i) = R(i).Data(1,1);
    M1(i) = R(i).Data(2,1);
    N1(i) = R(i).Data(3,1);
    L2(i) = R(i).Data(1,2);
    M2(i) = R(i).Data(2,2);
    N2(i) = R(i).Data(3,2);
    L3(i) = R(i).Data(1,3);
    M3(i) = R(i).Data(2,3);
    N3(i) = R(i).Data(3,3);
end

L1d = diff(L1)/300;
M1d = diff(M1)/300;
N1d = diff(N1)/300;
L2d = diff(L2)/300;
M2d = diff(M2)/300;
N2d = diff(N2)/300;
L3d = diff(L3)/300;
M3d = diff(M3)/300;
N3d = diff(N3)/300;

L1 = L1(2:end);
M1 = M1(2:end);
N1 = N1(2:end);
L2 = L2(2:end);
M2 = M2(2:end);
N2 = N2(2:end);
L3 = L3(2:end);
M3 = M3(2:end);
N3 = N3(2:end);

wx = L3.*L2d + M3.*M2d + N3.*N2d;
wy = L1.*L3d + M1.*M3d + N1.*N3d;
wz = L2.*L1d + M2.*M1d + N2.*N1d;

ax = diff(wx)/300;
ay = diff(wy)/300;
az = diff(wz)/300;


angVelo = [wx', wy', wz'];
angAccel = [ax', ay', az'];


end