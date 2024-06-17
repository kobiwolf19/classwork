clear
load("static0001.mat")
load("walk0002.mat")
addpath("/Users/kobi/Desktop/Biomechanical Methods")
fc = 10;
fs = 300;
walk = walk0002.Trajectories.Labeled.Data;
static = static0001.Trajectories.Labeled.Data;

walk(:,4,:) = [];

static(:,4,:) = [];

frames = walk0002.Frames;

time = linspace(0,frames/300, frames);

right_hip = [-150.9961 210.1181 947.1444];
left_hip = [-140.3969 377.7367 950.9822];

%% filter Data
filterstatic = struct('Data',{});
filterswalk = struct('Data',{});
for i = 1:height(static)
    xr_1 = squeeze(static(i,1,:));
    x = j211filter(xr_1,fc,fs);
    x = x/1000;
    yr_1 = squeeze(static(i,2,:));
    y = j211filter(yr_1,fc,fs);
    y = y/1000;
    zr_1 = squeeze(static(i,3,:));
    z = j211filter(zr_1,fc,fs);
    z = z/1000;
    fils = [x y z];
    filterstatic(i).Data = fils;
end
for i = 1:height(walk)
    xr_1 = squeeze(walk(i,1,:));
    x = j211filter(xr_1,fc,fs);
    x = x/1000;
    yr_1 = squeeze(walk(i,2,:));
    y = j211filter(yr_1,fc,fs);
    y = y/1000;
    zr_1 = squeeze(walk(i,3,:));
    z = j211filter(zr_1,fc,fs);
    z = z/1000;
    fils = [x y z];
    filterswalk(i).Data = fils;
end

%% 6.1

% right shank
r_knee = filterstatic(5).Data(2,:);
r_ankle = filterstatic(7).Data(2,:); %origin
r_shank = filterstatic(6).Data(2,:); %planar

right_shank = localCoordinateSystem(r_knee, r_ankle, r_shank);

% left shank
l_knee = filterstatic(11).Data(2,:);
l_ankle = filterstatic(13).Data(2,:); %origin
l_shank = filterstatic(12).Data(2,:); %planar

left_shank = localCoordinateSystem(l_knee, l_ankle, l_shank);

% right thigh

r_knee = filterstatic(5).Data(2,:);
r_thigh = filterstatic(4).Data(2,:); %planar

right_thigh = localCoordinateSystem(right_hip, r_knee, r_thigh);

% left thigh

l_knee = filterstatic(11).Data(2,:);
l_thigh = filterstatic(10).Data(2,:); %planar

left_thigh = localCoordinateSystem(left_hip, l_knee, l_thigh);

% right knee
p_o = filterstatic(16).Data(2,:) - filterstatic(5).Data(2,:);
pp = (right_thigh')*p_o';

%right ankle
p_o_ra = filterstatic(17).Data(2,:) - filterstatic(7).Data(2,:);
pp_ra = (right_shank')*p_o_ra';

%left knee
p_o_lk = filterstatic(18).Data(2,:) - filterstatic(11).Data(2,:);
pp_lk = (left_thigh')*p_o_lk';

%left ankle
p_o_la = filterstatic(19).Data(2,:) - filterstatic(13).Data(2,:);
pp_la = (left_shank')*p_o_la';

hipCoordinateSystem = struct('Data',{});

for i = 1:length(walk)

    %local coordinate system
    [r_hip_center(i,:), l_hip_center(i,:), hipCoordinateSystem(i).Data] = hipCenter(filterswalk(1).Data(i,:),filterswalk(3).Data(i,:),filterswalk(2).Data(i,:));

    %right shank
    r_knee = filterswalk(5).Data(i,:);
    r_ankle = filterswalk(7).Data(i,:); %origin
    r_shank = filterswalk(6).Data(i,:); %planar

    right_shank_walk = localCoordinateSystem(r_knee, r_ankle, r_shank);

    % left shank
    l_knee = filterswalk(11).Data(i,:);
    l_ankle = filterswalk(13).Data(i,:); %origin
    l_shank = filterswalk(12).Data(i,:); %planar

    left_shank_walk = localCoordinateSystem(l_knee, l_ankle, l_shank);

    % right thigh
    r_knee = filterswalk(5).Data(i,:);
    r_thigh = filterswalk(4).Data(i,:); %planar

    right_thigh_walk = localCoordinateSystem(r_hip_center(i,:), r_knee, r_thigh);

    % left thigh
    l_knee = filterswalk(11).Data(i,:);
    l_thigh = filterswalk(10).Data(i,:); %planar

    left_thigh_walk = localCoordinateSystem(l_hip_center(i,:), l_knee, l_thigh);


    % P prime
    %right knee
    p_o = filterswalk(16).Data(i,:) - filterswalk(5).Data(i,:);
    pp = (right_thigh')*p_o';

    %right ankle
    p_o_ra = filterswalk(17).Data(i,:) - filterswalk(7).Data(i,:);
    pp_ra = (right_shank')*p_o_ra';

    %left knee
    p_o_lk = filterswalk(18).Data(i,:) - filterswalk(11).Data(i,:);
    pp_lk = (left_thigh')*p_o_lk';

    %left ankle
    p_o_la = filterswalk(19).Data(i,:) - filterswalk(13).Data(i,:);
    pp_la = (left_shank')*p_o_la';
    
    % final equation
    P_rightthigh = right_thigh_walk*pp + r_knee';
    P_r_mknee_final(i,:) = P_rightthigh';

    P_leftthigh = left_thigh_walk*pp_lk + l_knee';
    P_l_mknee_final(i,:) = P_leftthigh';

    P_rightshank = right_shank*pp_ra + r_ankle';
    P_r_mankle_final(i,:) = P_rightshank';

    P_leftshank = left_shank*pp_la + l_ankle';
    P_l_mankle_final(i,:) = P_leftshank';

end


l_knee_center = zeros(847,3);
r_knee_center = zeros(847,3);
l_ankle_center = zeros(847,3);
r_ankle_center = zeros(847,3);



for i = 1:length(walk)
    l_knee_center(i,1:3) = (filterswalk(11).Data(i,:) + P_l_mknee_final(i,:))/2;
    r_knee_center(i,1:3) = (filterswalk(5).Data(i,:) + P_r_mknee_final(i,:))/2;
    
    l_ankle_center(i,1:3) = (filterswalk(13).Data(i,:) + P_l_mankle_final(i,:))/2;
    r_ankle_center(i,1:3) = (filterswalk(7).Data(i,:) + P_r_mankle_final(i,:))/2;
end


plot(l_knee_center)
legend('x','y','z')
grid on
title('Left Knee Center')

plot(r_knee_center)
legend('x','y','z')
grid on
title('Right Knee Center')

plot(l_ankle_center)
legend('x','y','z')
grid on
title('Left Ankle Center')

plot(r_ankle_center)
legend('x','y','z')
grid on
title('Right Ankle Center')


%% 6.2

% Rotation Matrix
anatomicalCoor_lthigh = struct('Data',{});
anatomicalCoor_rthigh = struct('Data',{});
anatomicalCoor_lshank = struct('Data',{});
anatomicalCoor_rshank = struct('Data',{});
anatomicalCoor_lfoot = struct('Data',{});
anatomicalCoor_rfoot = struct('Data',{});

for i = 1:length(walk)

    % Left Thigh
    lthigh_length = norm(l_knee_center(i,:) - l_hip_center(i,:));
    IS_axis_lthigh = ((l_knee_center(i,:) - l_hip_center(i,:))/lthigh_length)';

    temp_length = norm(l_knee_center(i,:) - filterswalk(11).Data(i,:));
    temp_lthigh = ((l_knee_center(i,:) - filterswalk(11).Data(i,:))/temp_length)';

    AP_axis_lthigh = cross(temp_lthigh, IS_axis_lthigh);

    ML_axis_lthigh = cross(IS_axis_lthigh, AP_axis_lthigh);

    anatomicalCoor_lthigh(i).Data = horzcat(IS_axis_lthigh, AP_axis_lthigh, ML_axis_lthigh);

    % Right Thigh
    rthigh_length = norm(r_knee_center(i,:) - r_hip_center(i,:));
    IS_axis_rthigh = ((r_knee_center(i,:) - r_hip_center(i,:))/rthigh_length)';

    temp_length = norm(filterswalk(5).Data(i,:) - r_knee_center(i,:));
    temp_rthigh = ((filterswalk(5).Data(i,:) - r_knee_center(i,:))/temp_length)';

    AP_axis_rthigh = cross(temp_rthigh, IS_axis_rthigh);

    ML_axis_rthigh = cross(IS_axis_rthigh, AP_axis_rthigh);

    anatomicalCoor_rthigh(i).Data = horzcat(IS_axis_rthigh, AP_axis_rthigh, ML_axis_rthigh);

    % Left shank
    lshank_length = norm(l_ankle_center(i,:) - l_knee_center(i,:));
    IS_axis_lshank = ((l_ankle_center(i,:) - l_knee_center(i,:))/lshank_length)';

    temp_length = norm(l_ankle_center(i,:) - filterswalk(13).Data(i,:));
    temp_lshank = ((l_ankle_center(i,:) - filterswalk(13).Data(i,:))/temp_length)';

    AP_axis_lshank = cross(temp_lshank, IS_axis_lshank);

    ML_axis_lshank = cross(IS_axis_lshank, AP_axis_lshank);

    anatomicalCoor_lshank(i).Data = horzcat(IS_axis_lshank, AP_axis_lshank, ML_axis_lshank);

    % Right shank
    rshank_length = norm(r_ankle_center(i,:) - r_knee_center(i,:));
    IS_axis_rshank = ((r_ankle_center(i,:) - r_knee_center(i,:))/rshank_length)';

    temp_length = norm(filterswalk(7).Data(i,:) - r_ankle_center(i,:));
    temp_rshank = ((filterswalk(7).Data(i,:) - r_ankle_center(i,:))/temp_length)';

    AP_axis_rshank = cross(temp_rshank, IS_axis_rshank);

    ML_axis_rshank = cross(IS_axis_rshank, AP_axis_rshank);

    anatomicalCoor_rshank(i).Data = horzcat(IS_axis_rshank, AP_axis_rshank, ML_axis_rshank);

    % Left foot
    left_foot_length = norm(filterswalk(15).Data(i,:) - filterswalk(14).Data(i,:));
    left_foot_Y = ((filterswalk(15).Data(i,:) - filterswalk(14).Data(i,:))/left_foot_length)';

    temp_length = norm(l_ankle_center(i,:) - filterswalk(13).Data(i,:));
    temp_lfoot = ((l_ankle_center(i,:) - filterswalk(13).Data(i,:))/temp_length)';

    left_foot_X = cross(left_foot_Y, temp_lfoot);

    left_foot_Z = cross(left_foot_X, left_foot_Y);

    anatomicalCoor_lfoot(i).Data = horzcat(left_foot_X, left_foot_Y, left_foot_Z);

    % Right Foot
    right_foot_length = norm(filterswalk(9).Data(i,:) - filterswalk(8).Data(i,:));
    right_foot_Y = ((filterswalk(9).Data(i,:) - filterswalk(8).Data(i,:))/right_foot_length)';

    temp_length = norm(filterswalk(7).Data(i,:) - r_ankle_center(i,:));
    temp_rfoot = ((filterswalk(7).Data(i,:) - r_ankle_center(i,:))/temp_length)';

    right_foot_X = cross(right_foot_Y, temp_rfoot);

    right_foot_Z = cross(right_foot_X, right_foot_Y);

    anatomicalCoor_rfoot(i).Data = horzcat(right_foot_X, right_foot_Y, right_foot_Z);
end

%% Assignment 7

JointAngle_rHip = zeros(847,3);
JointAngle_lHip = zeros(847,3);
JointAngle_rKnee = zeros(847,3);
JointAngle_lKnee = zeros(847,3);
JointAngle_rAnkle = zeros(847,3);
JointAngle_lAnkle = zeros(847,3);

for i = 1:frames

    [a_rHip, b_rHip, y_rHip] = jointAngles(hipCoordinateSystem(i).Data, anatomicalCoor_rthigh(i).Data);
    JointAngle_rHip(i,1:3) = [a_rHip, b_rHip, y_rHip];

    [a_lHip, b_lHip, y_lHip] = jointAngles(hipCoordinateSystem(i).Data, anatomicalCoor_lthigh(i).Data);
    JointAngle_lHip(i,1:3) = [a_lHip, b_lHip, y_lHip];

    [a_rShank, b_rShank, y_rShank] = jointAngles(anatomicalCoor_rthigh(i).Data, anatomicalCoor_rshank(i).Data);
    JointAngle_rKnee(i,1:3) = [a_rShank, b_rShank, y_rShank];

    [a_lShank, b_lShank, y_lShank] = jointAngles(anatomicalCoor_lthigh(i).Data, anatomicalCoor_lshank(i).Data);
    JointAngle_lKnee(i,1:3) = [a_lShank, b_lShank, y_lShank];

    [a_rAnkle, b_rAnkle, y_rAnkle] = jointAngles(anatomicalCoor_rshank(i).Data, anatomicalCoor_rfoot(i).Data);
    JointAngle_rAnkle(i,1:3) = [a_rAnkle, b_rAnkle, y_rAnkle];

    [a_lAnkle, b_lAnkle, y_lAnkle] = jointAngles(anatomicalCoor_lshank(i).Data, anatomicalCoor_lfoot(i).Data);
    JointAngle_lAnkle(i,1:3) = [a_lAnkle, b_lAnkle, y_lAnkle];

end



plot(JointAngle_rHip(:,1))
hold on
plot(JointAngle_rHip(:,2))
plot(JointAngle_rHip(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Right Hip Euler Angles')
grid on
hold off

plot(JointAngle_lHip(:,1))
hold on
plot(JointAngle_lHip(:,2))
plot(JointAngle_lHip(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Left Hip Euler Angles')
grid on
hold off

plot(JointAngle_rKnee(:,1))
hold on
plot(JointAngle_rKnee(:,2))
plot(JointAngle_rKnee(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Right Knee Euler Angles')
grid on
hold off

plot(JointAngle_lKnee(:,1))
hold on
plot(JointAngle_lKnee(:,2))
plot(JointAngle_lKnee(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Left Knee Euler Angles')
grid on
hold off

plot(JointAngle_rAnkle(:,1))
hold on
plot(JointAngle_rAnkle(:,2))
plot(JointAngle_rAnkle(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Right Ankle Euler Angles')
grid on
hold off

plot(JointAngle_lAnkle(:,1))
hold on
plot(JointAngle_lAnkle(:,2))
plot(JointAngle_lAnkle(:,3))
legend('flexion/extension','ab/adduction','internal/external rotation')
title('Left Ankle Euler Angles')
grid on
hold off

%% assignment 8

Hel_rHip = zeros(847,3);
Hel_lHip = zeros(847,3);
Hel_rKnee = zeros(847,3);
Hel_lKnee = zeros(847,3);
Hel_rAnkle = zeros(847,3);
Hel_lAnkle = zeros(847,3);

for i = 1:frames

    K = HelicalAngles(hipCoordinateSystem(i).Data, anatomicalCoor_rthigh(i).Data);
    Hel_rHip(i,1:3) = [K(1), K(2), K(3)];

    K = HelicalAngles(hipCoordinateSystem(i).Data, anatomicalCoor_lthigh(i).Data);
    Hel_lHip(i,1:3) = [K(1), K(2), K(3)];

    K = HelicalAngles(anatomicalCoor_rthigh(i).Data, anatomicalCoor_rshank(i).Data);
    Hel_rKnee(i,1:3) = [K(1), K(2), K(3)];

    K = HelicalAngles(anatomicalCoor_lthigh(i).Data, anatomicalCoor_lshank(i).Data);
    Hel_lKnee(i,1:3) = [K(1), K(2), K(3)];

    K = HelicalAngles(anatomicalCoor_rshank(i).Data, anatomicalCoor_rfoot(i).Data);
    Hel_rAnkle(i,1:3) = [K(1), K(2), K(3)];

    K = HelicalAngles(anatomicalCoor_lshank(i).Data, anatomicalCoor_lfoot(i).Data);
    Hel_lAnkle(i,1:3) = [K(1), K(2), K(3)];

end
   

plot(Hel_rHip(:,1))
hold on
plot(Hel_rHip(:,2))
plot(Hel_rHip(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Right Hip Helical Angles')
grid on
hold off

plot(Hel_lHip(:,1))
hold on
plot(Hel_lHip(:,2))
plot(Hel_lHip(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Left Hip Helical Angles')
grid on
hold off

plot(Hel_rKnee(:,1))
hold on
plot(Hel_rKnee(:,2))
plot(Hel_rKnee(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Right Knee Helical Angles')
grid on
hold off

plot(Hel_lKnee(:,1))
hold on
plot(Hel_lKnee(:,2))
plot(Hel_lKnee(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Left Knee Helical Angles')
grid on
hold off

plot(Hel_rAnkle(:,1))
hold on
plot(Hel_rAnkle(:,2))
plot(Hel_rAnkle(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Right Ankle Helical Angles')
grid on
hold off

plot(Hel_lAnkle(:,1))
hold on
plot(Hel_lAnkle(:,2))
plot(Hel_lAnkle(:,3))
legend('internal/external rotation', 'ab/adduction', 'flexion/extension')
title('Left Ankle Helical Angles')
grid on
hold off


%% Assignment 9

% 9.1

COM_rThigh = (r_knee_center - r_hip_center)*(.433) + r_hip_center;
COM_lThigh = (l_knee_center - l_hip_center)*(.433) + l_hip_center;

COM_rShank = (r_ankle_center - r_knee_center) * 0.433 + r_knee_center;
COM_lShank = (l_ankle_center - l_knee_center) * 0.433 + l_knee_center;

COM_rFoot = (filterswalk(9).Data - r_ankle_center) * 0.5 + r_ankle_center;
COM_lFoot = (filterswalk(15).Data - l_ankle_center) * 0.5 + l_ankle_center;

% 9.2 vectors
rF_TP = (r_hip_center - COM_rThigh);
lF_TP = (l_hip_center - COM_lThigh);

rF_TD = (r_knee_center - COM_rThigh);
lF_TD = (l_knee_center - COM_lThigh);

rF_SP = (r_knee_center - COM_rShank);
lF_SP = (l_knee_center - COM_lShank);

rF_SD = (r_ankle_center - COM_rShank);
lF_SD = (l_ankle_center - COM_lShank);

rF_FP = (r_ankle_center - COM_rFoot);
lF_FP = (l_ankle_center - COM_lFoot);


COP_R = walk0002.Force(2).COP';
COP_R = downsample(COP_R,4);
COP_R = COP_R/1000;

COP_L = walk0002.Force(3).COP'; 
COP_L = downsample(COP_L,4);
COP_L = COP_L/1000;

index_R = find(COP_R, 1, 'first');
index_L = find(COP_L, 1, 'first');

F_D_R = (COP_R - COM_rFoot);
F_D_L = (COP_L - COM_lFoot);

plot(rF_FP)

%% Assignment 10

vel_COM_RF = diff(COM_rFoot)/300;
acc_COM_RF = diff(vel_COM_RF)/300;

vel_COM_LF = diff(COM_lFoot)/300;
acc_COM_LF = diff(vel_COM_LF)/300;

veL_COM_RS = diff(COM_rShank)/300;
acc_COM_RS = diff(veL_COM_RS)/300;

veL_COM_LS = diff(COM_lShank)/300;
acc_COM_LS = diff(veL_COM_LS)/300;

vel_COM_RT = diff(COM_rThigh)/300;
acc_COM_RT = diff(vel_COM_RT)/300;

vel_COM_LT = diff(COM_lThigh)/300;
acc_COM_LT = diff(vel_COM_LT)/300;


% ground reaction forces
GRF_R = walk0002.Force(2).Force';
GRF_R = downsample(GRF_R,4);
GRF_R(1:2,:) = [];

GRF_L = walk0002.Force(3).Force';
GRF_L = downsample(GRF_L,4);
GRF_L(1:2,:) = [];

% distal foot force
rFoot_Distal = zeros(length(GRF_R),3);
rFoot_Distal(:,1) = -GRF_R(:,1);
rFoot_Distal(:,2) = -GRF_R(:,2);
rFoot_Distal(:,3) = -GRF_R(:,3);

lFoot_Force = zeros(length(GRF_L),3);
lFoot_Force(:,1) = -GRF_L(:,1);
lFoot_Force(:,2) = -GRF_L(:,2);
lFoot_Force(:,3) = -GRF_L(:,3);

% proximal foot force
foot_mass = 61.2*0.0137;
Rfoot_prox = zeros(length(GRF_R),3);
Rfoot_prox(:,1) = foot_mass*acc_COM_RF(:,1) - rFoot_Distal(:,1);
Rfoot_prox(:,2) = foot_mass*acc_COM_RF(:,2) - rFoot_Distal(:,2);
Rfoot_prox(:,3) = foot_mass*acc_COM_RF(:,3) - rFoot_Distal(:,3) - foot_mass*9.81;

Lfoot_prox = zeros(length(GRF_R),3);
Lfoot_prox(:,1) = foot_mass*acc_COM_LF(:,1) - lFoot_Force(:,1);
Lfoot_prox(:,2) = foot_mass*acc_COM_LF(:,2) - lFoot_Force(:,2);
Lfoot_prox(:,3) = foot_mass*acc_COM_LF(:,3) - lFoot_Force(:,3) - foot_mass*9.81;

% distal shank force
Rshank_dist = -Rfoot_prox;
Lshank_dist = -Lfoot_prox;

% proximal shank force
shank_mass = 61.2*0.0433;
Rshank_prox = zeros(length(GRF_R),3);
Rshank_prox(:,1) = shank_mass*acc_COM_RS(:,1) - Rshank_dist(:,1);
Rshank_prox(:,2) = shank_mass*acc_COM_RS(:,2) - Rshank_dist(:,2);
Rshank_prox(:,3) = shank_mass*acc_COM_RS(:,3) - Rshank_dist(:,3) - shank_mass*9.81;

Lshank_prox = zeros(length(GRF_R),3);
Lshank_prox(:,1) = shank_mass*acc_COM_LS(:,1) - Lshank_dist(:,1);
Lshank_prox(:,2) = shank_mass*acc_COM_LS(:,2) - Lshank_dist(:,2);
Lshank_prox(:,3) = shank_mass*acc_COM_LS(:,3) - Lshank_dist(:,3) - shank_mass*9.81;

% distal thigh force
Rthigh_dist = -Rshank_prox;
Lthigh_dist = -Lshank_prox;

% proximal thigh force
thigh_mass = 61.2*0.14165;
Rthigh_prox = zeros(length(GRF_R),3);
Rthigh_prox(:,1) = thigh_mass*acc_COM_RT(:,1) - Rthigh_dist(:,1);
Rthigh_prox(:,2) = thigh_mass*acc_COM_RT(:,2) - Rthigh_dist(:,2);
Rthigh_prox(:,3) = thigh_mass*acc_COM_RT(:,3) - Rthigh_dist(:,3) - thigh_mass*9.81;

Lthigh_prox = zeros(length(GRF_R),3);
Lthigh_prox(:,1) = thigh_mass*acc_COM_LT(:,1) - Lthigh_dist(:,1);
Lthigh_prox(:,2) = thigh_mass*acc_COM_LT(:,2) - Lthigh_dist(:,2);
Lthigh_prox(:,3) = thigh_mass*acc_COM_LT(:,3) - Lthigh_dist(:,3) - thigh_mass*9.81;


plot(Rthigh_prox)
grid on
legend('x','y','z')
title('Right Thigh Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

plot(Lthigh_prox)
grid on
legend('x','y','z')
title('Left Thigh Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

plot(Rshank_prox)
grid on
legend('x','y','z')
title('Right Shank Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

plot(Lshank_prox)
grid on
legend('x','y','z')
title('Left Shank Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

plot(Rfoot_prox)
grid on
legend('x','y','z')
title('Right Foot Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

plot(Lfoot_prox)
grid on
legend('x','y','z')
title('Left Foot Proximal Forces')
ylabel('Position (meters)')
xlabel('Frames')

%% Assignment 11
% 11.1

global_Rfoot_MCG = cross(rF_FP(3:847,:),Rfoot_prox) + cross(F_D_R(3:847,:),rFoot_Distal);
global_Lfoot_MCG = cross(lF_FP(3:847,:),Lfoot_prox) + cross(F_D_L(3:847,:),lFoot_Force);

global_Rshank_MCG = cross(rF_SP(3:847,:),(Rshank_prox)) + cross(rF_SD(3:847,:),Rshank_dist);
global_Lshank_MCG = cross(lF_SP(3:847,:),(Lshank_prox)) + cross(lF_SD(3:847,:),Lshank_dist);

global_Rthigh_MCG = cross(rF_TP(3:847,:),(Rthigh_prox)) + cross(rF_TD(3:847,:),Rthigh_dist);
global_Lthigh_MCG = cross(lF_TP(3:847,:),(Lthigh_prox)) + cross(lF_TD(3:847,:),Lthigh_dist);

plot(global_Rfoot_MCG)
grid on
legend('x','y','z')
title('Right Foot Global Moments')
ylabel('Moment')
xlabel('Frames')

plot(global_Lfoot_MCG)
grid on
legend('x','y','z')
title('Left Foot Global Moments')
ylabel('Moment')
xlabel('Frames')

plot(global_Rshank_MCG)
grid on
legend('x','y','z')
title('Right Shank Global Moments')
ylabel('Moment')
xlabel('Frames')

plot(global_Lshank_MCG)
grid on
legend('x','y','z')
title('Left Shank Global Moments')
ylabel('Moment')
xlabel('Frames')

plot(global_Rthigh_MCG)
grid on
legend('x','y','z')
title('Right Thigh Global Moments')
ylabel('Moment')
xlabel('Frames')

plot(global_Lthigh_MCG)
grid on
legend('x','y','z')
title('Left Thigh Global Moments')
ylabel('Moment')
xlabel('Frames')

global_Rfoot_MCG(index_R,:)
global_Lfoot_MCG(index_L,:)
global_Rshank_MCG(index_R,:)
global_Lshank_MCG(index_L,:)
global_Rthigh_MCG(index_R,:)
global_Lthigh_MCG(index_L,:)


% 11.2
inertia_foot(1) = 0;
inertia_foot(2) = 0.0061*((61.2*1.753^2)/(73.59*1.755^2));
inertia_foot(3) = 0.0061*((61.2*1.753^2)/(73.59*1.755^2));

inertia_shank(1) = 0.0001*((61.2*1.753^2)/(73.59*1.755^2));
inertia_shank(2) = 0.0414*((61.2*1.753^2)/(73.59*1.755^2));
inertia_shank(3) = 0.0414*((61.2*1.753^2)/(73.59*1.755^2));

inertia_thigh(1) = 0.0005*((61.2*1.753^2)/(73.59*1.755^2));
inertia_thigh(2) = 0.0878*((61.2*1.753^2)/(73.59*1.755^2));
inertia_thigh(3) = 0.0878*((61.2*1.753^2)/(73.59*1.755^2));

%% 12.1

Lthigh_angvelo_temp = zeros(845,3);
Rthigh_angvelo_temp = zeros(845,3);
Lshank_angvelo_temp = zeros(845,3);
Rshank_angvelo_temp = zeros(845,3);
Lfoot_angvelo_temp = zeros(845,3);
Rfoot_angvelo_temp = zeros(845,3);

for i = 1:frames - 2
    Lthigh_angvelo_temp(i,1:3) = anatomicalCoor_lthigh(i).Data * global_Lthigh_MCG(i,1:3)';
    Rthigh_angvelo_temp(i,1:3) = anatomicalCoor_rthigh(i).Data * global_Rthigh_MCG(i,1:3)';
    Lshank_angvelo_temp(i,1:3) = anatomicalCoor_lshank(i).Data * global_Lshank_MCG(i,1:3)';
    Rshank_angvelo_temp(i,1:3) = anatomicalCoor_rshank(i).Data * global_Rshank_MCG(i,1:3)';
    Lfoot_angvelo_temp(i,1:3) = anatomicalCoor_lfoot(i).Data * global_Lfoot_MCG(i,1:3)';
    Rfoot_angvelo_temp(i,1:3) = anatomicalCoor_rfoot(i).Data * global_Rfoot_MCG(i,1:3)';
end

figure(1)
subplot(3,2,1)
plot(Lthigh_angvelo_temp)
grid on
legend('x','y','z')
title('Left Thigh Local Moments')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,2)
plot(Rthigh_angvelo_temp)
grid on
legend('x','y','z')
title('Right Thigh Local Moments')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,3)
plot(Lshank_angvelo_temp)
grid on
legend('x','y','z')
title('Left Shank Local Moments')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,4)
plot(Rshank_angvelo_temp)
grid on
legend('x','y','z')
title('Right Shank Local Moments')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,5)
plot(Lfoot_angvelo_temp)
grid on
legend('x','y','z')
title('Left Foot Local Moments')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,6)
plot(Rfoot_angvelo_temp)
grid on
legend('x','y','z')
title('Right Foot Local Moments')
ylabel('Moment')
xlabel('Frames')

%% 12.2
Moment_R = walk0002.Force(2).Moment';
Moment_R = downsample(Moment_R,4);
Moment_R = Moment_R/1000;
Moment_L = walk0002.Force(3).Moment';
Moment_L = downsample(Moment_L,4);
Moment_L = Moment_L/1000;

Lfoot_moment_FP = zeros(845,3);
Rfoot_moment_FP = zeros(845,3);

% finding moments with the force plates
for i = 1:frames - 2
    Lfoot_moment_FP(i,1:3) = anatomicalCoor_lfoot(i).Data * Moment_L(i,1:3)';
    Rfoot_moment_FP(i,1:3) = anatomicalCoor_rfoot(i).Data * Moment_R(i,1:3)';
end

figure(2)
subplot(2,1,1)
plot(Lfoot_moment_FP)
grid on
legend('x','y','z')
title('Left Foot Local Moments from Force Plates')
ylabel('Moment')
xlabel('Frames')

subplot(2,1,2)
plot(Rfoot_moment_FP)
grid on
legend('x','y','z')
title('Right Foot Local Moments from Force Plates')
ylabel('Moment')
xlabel('Frames')

%% 12.3 angular velocities and accelerations

[Rthigh_angVelo, Rthigh_angAccel] = getAngVelo(anatomicalCoor_rthigh);
[Lthigh_angVelo, Lthigh_angAccel] = getAngVelo(anatomicalCoor_lthigh);
[Rshank_angVelo, Rshank_angAccel] = getAngVelo(anatomicalCoor_rshank);
[Lshank_angVelo, Lshank_angAccel] = getAngVelo(anatomicalCoor_lshank);
[Rfoot_angVelo, Rfoot_angAccel] = getAngVelo(anatomicalCoor_rfoot);
[Lfoot_angVelo, Lfoot_angAccel] = getAngVelo(anatomicalCoor_lfoot);

figure(3)
subplot(3,2,1)
plot(Rthigh_angVelo)
legend('x','y','z')
grid on
title('Right Thigh Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,2)
plot(Rthigh_angAccel)
legend('x','y','z')
grid on
title('Right Thigh Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')

subplot(3,2,3)
plot(Lthigh_angVelo)
legend('x','y','z')
grid on
title('Left Thigh Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,4)
plot(Lthigh_angAccel)
legend('x','y','z')
grid on
title('Left Thigh Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')

subplot(3,2,5)
plot(Rshank_angVelo)
legend('x','y','z')
grid on
title('Right Shank Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,6)
plot(Rshank_angAccel)
legend('x','y','z')
grid on
title('Right Shank Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')

figure(6)
subplot(3,2,1)
plot(Lshank_angVelo)
legend('x','y','z')
grid on
title('Left Shank Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,2)
plot(Lshank_angAccel)
legend('x','y','z')
grid on
title('Left Shank Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')

subplot(3,2,3)
plot(Rfoot_angVelo)
legend('x','y','z')
grid on
title('Right Foot Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,4)
plot(Rfoot_angAccel)
legend('x','y','z')
grid on
title('Right Foot Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')

subplot(3,2,5)
plot(Lfoot_angVelo)
legend('x','y','z')
grid on
title('Left Foot Angular Velocity')
ylabel('Angular Velocity')
xlabel('Frames')

subplot(3,2,6)
plot(Lfoot_angAccel)
legend('x','y','z')
grid on
title('Left Foot Angular Acceleration')
ylabel('Angular Acceleration')
xlabel('Frames')


%% Assignment 13

% 13.1 proximal segment moments

Mp_RFoot(:,1) = -Moment_R(3:end,1) - global_Rfoot_MCG(:,1) + inertia_foot(1) .* Rfoot_angAccel(:,1) - Rfoot_angVelo(2:end,2) .* Rfoot_angVelo(2:end,3) .* (inertia_foot(2) - inertia_foot(3));
Mp_RFoot(:,2) = -Moment_R(3:end,2) - global_Rfoot_MCG(:,2) + inertia_foot(2) .* Rfoot_angAccel(:,2) - Rfoot_angVelo(2:end,1) .* Rfoot_angVelo(2:end,3) .* (inertia_foot(3) - inertia_foot(1));
Mp_RFoot(:,3) = -Moment_R(3:end,3) - global_Rfoot_MCG(:,3) + inertia_foot(3) .* Rfoot_angAccel(:,3) - Rfoot_angVelo(2:end,2) .* Rfoot_angVelo(2:end,1) .* (inertia_foot(1) - inertia_foot(2));

Mp_LFoot(:,1) = -Moment_L(3:end,1) - global_Lfoot_MCG(:,1) + inertia_foot(1) .* Lfoot_angAccel(:,1) - Lfoot_angVelo(2:end,2) .* Lfoot_angVelo(2:end,3) .* (inertia_foot(2) - inertia_foot(3));
Mp_LFoot(:,2) = -Moment_L(3:end,2) - global_Lfoot_MCG(:,2) + inertia_foot(2) .* Lfoot_angAccel(:,2) - Lfoot_angVelo(2:end,1) .* Lfoot_angVelo(2:end,3) .* (inertia_foot(3) - inertia_foot(1));
Mp_LFoot(:,3) = -Moment_L(3:end,3) - global_Lfoot_MCG(:,3) + inertia_foot(3) .* Lfoot_angAccel(:,3) - Lfoot_angVelo(2:end,2) .* Lfoot_angVelo(2:end,1) .* (inertia_foot(1) - inertia_foot(2));


Mp_RShank(:,1) = -Mp_RFoot(:,1) - global_Rshank_MCG(:,1) + inertia_shank(1) .* Rshank_angAccel(:,1) - Rshank_angVelo(2:end,2) .* Rshank_angVelo(2:end,3) .* (inertia_shank(2) - inertia_shank(3));
Mp_RShank(:,2) = -Mp_RFoot(:,2) - global_Rshank_MCG(:,2) + inertia_shank(2) .* Rshank_angAccel(:,2) - Rshank_angVelo(2:end,1) .* Rshank_angVelo(2:end,3) .* (inertia_shank(3) - inertia_shank(1));
Mp_RShank(:,3) = -Mp_RFoot(:,3) - global_Rshank_MCG(:,3) + inertia_shank(3) .* Rshank_angAccel(:,3) - Rshank_angVelo(2:end,2) .* Rshank_angVelo(2:end,1) .* (inertia_shank(1) - inertia_shank(2));

Mp_LShank(:,1) = -Mp_LFoot(:,1) - global_Lshank_MCG(:,1) + inertia_shank(1) .* Lshank_angAccel(:,1) - Lshank_angVelo(2:end,2) .* Lshank_angVelo(2:end,3) .* (inertia_shank(2) - inertia_shank(3));
Mp_LShank(:,2) = -Mp_LFoot(:,2) - global_Lshank_MCG(:,2) + inertia_shank(2) .* Lshank_angAccel(:,2) - Lshank_angVelo(2:end,1) .* Lshank_angVelo(2:end,3) .* (inertia_shank(3) - inertia_shank(1));
Mp_LShank(:,3) = -Mp_LFoot(:,3) - global_Lshank_MCG(:,3) + inertia_shank(3) .* Lshank_angAccel(:,3) - Lshank_angVelo(2:end,2) .* Lshank_angVelo(2:end,1) .* (inertia_shank(1) - inertia_shank(2));


Mp_RThigh(:,1) = -Mp_RShank(:,1) - global_Rthigh_MCG(:,1) + inertia_thigh(1) .* Rthigh_angAccel(:,1) - Rthigh_angVelo(2:end,2) .* Rthigh_angVelo(2:end,3) .* (inertia_thigh(2) - inertia_thigh(3));
Mp_RThigh(:,2) = -Mp_RShank(:,2) - global_Rthigh_MCG(:,2) + inertia_thigh(2) .* Rthigh_angAccel(:,2) - Rthigh_angVelo(2:end,1) .* Rthigh_angVelo(2:end,3) .* (inertia_thigh(3) - inertia_thigh(1));
Mp_RThigh(:,3) = -Mp_RShank(:,3) - global_Rthigh_MCG(:,3) + inertia_thigh(3) .* Rthigh_angAccel(:,3) - Rthigh_angVelo(2:end,2) .* Rthigh_angVelo(2:end,1) .* (inertia_thigh(1) - inertia_thigh(2));

Mp_LThigh(:,1) = -Mp_LShank(:,1) - global_Lthigh_MCG(:,1) + inertia_thigh(1) .* Lthigh_angAccel(:,1) - Lthigh_angVelo(2:end,2) .* Lthigh_angVelo(2:end,3) .* (inertia_thigh(2) - inertia_thigh(3));
Mp_LThigh(:,2) = -Mp_LShank(:,2) - global_Lthigh_MCG(:,2) + inertia_thigh(2) .* Lthigh_angAccel(:,2) - Lthigh_angVelo(2:end,1) .* Lthigh_angVelo(2:end,3) .* (inertia_thigh(3) - inertia_thigh(1));
Mp_LThigh(:,3) = -Mp_LShank(:,3) - global_Lthigh_MCG(:,3) + inertia_thigh(3) .* Lthigh_angAccel(:,3) - Lthigh_angVelo(2:end,2) .* Lthigh_angVelo(2:end,1) .* (inertia_thigh(1) - inertia_thigh(2));

figure(4)
subplot(3,2,1)
plot(Mp_RFoot)
legend('x','y','z')
grid on
grid on
title('Right Foot Proximal Moment')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,2)
plot(Mp_LFoot)
legend('x','y','z')
grid on
grid on
title('Left Foot Proximal Moment')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,3)
plot(Mp_RShank)
legend('x','y','z')
grid on
grid on
title('Right Shank Proximal Moment')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,4)
plot(Mp_LShank)
legend('x','y','z')
grid on
grid on
title('Left Shank Proximal Moment')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,5)
plot(Mp_RThigh)
legend('x','y','z')
grid on
grid on
title('Right Thigh Proximal Moment')
ylabel('Moment')
xlabel('Frames')

subplot(3,2,6)
plot(Mp_LThigh)
legend('x','y','z')
grid on
grid on
title('Left Thigh Proximal Moment')
ylabel('Moment')
xlabel('Frames')


%% 13.2 hip, knee, and ankle powers in the sagittal plane

power_Rfoot = Mp_RFoot .* Rfoot_angVelo(2:end,:);
power_Lfoot = Mp_LFoot .* Lfoot_angVelo(2:end,:);
power_Rshank = Mp_RShank .* Rshank_angVelo(2:end,:);
power_Lshank = Mp_LShank .* Lshank_angVelo(2:end,:);
power_Rthigh = Mp_RThigh .* Rthigh_angVelo(2:end,:);
power_Lthigh = Mp_LThigh .* Lthigh_angVelo(2:end,:);


figure(5)
subplot(3,2,1)
plot(power_Rfoot(:,1))
grid on
title('Right Foot Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

subplot(3,2,2)
plot(power_Lfoot(:,1))
grid on
title('Left Foot Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

subplot(3,2,3)
plot(power_Rshank(:,1))
grid on
title('Right Shank Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

subplot(3,2,4)
plot(power_Lshank(:,1))
grid on
title('Left Shank Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

subplot(3,2,5)
plot(power_Rthigh(:,1))
grid on
title('Right Thigh Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

subplot(3,2,6)
plot(power_Lthigh(:,1))
grid on
title('Left Thigh Power in the Sagittal Plane')
ylabel('Power')
xlabel('Frames')

