%% head calculations
VW_L_head = normcdf((log(169)-7.45231)/0.73998);
VW_R_head = normcdf((log(118)-7.45231)/0.73998);

%% 1b driver neck
load v11660_009.txt
L_Z_force = v11660_009;
load v11660_011.txt
L_Y_moment = v11660_011;

L_Z_force = j211filtfilt(600,2900,L_Z_force(:,2));
L_Y_moment = j211filtfilt(600,2900,L_Y_moment(:,2));

L_Z_force = L_Z_force/1000;

NIJ = zeros(length(L_Y_moment),1);

check_maxNIJ = 0;
NIJ_ind_L = 0;
load_case = '';

for i = 1:length(L_Z_force)
    Fz = L_Z_force(i);
    My = L_Y_moment(i);

    if Fz > 0
        CVForce = 6806/1000; %tension
    else
        CVForce = 6160/1000; %compression
    end
    if My > 0
        CVMoment = 310; %flexion
    else
        CVMoment = 135; %extension
    end

    NIJ(i) = abs(Fz/CVForce) + abs(My/CVMoment);

    if check_maxNIJ < NIJ(i)
        check_maxNIJ = NIJ(i);
        NIJ_ind_L = i;
        if CVForce == 6.806 && CVMoment == 310
            load_case = 'tension/flexion';
        elseif CVForce == 6.16 && CVMoment == 310
            load_case = 'compression/flexion';
        elseif CVForce == 6.806 && CVMoment == 135
            load_case = 'tension/extension';
        else
            load_case = 'compression/extension';
        end
    end
end

maxNIJ_L = max(NIJ);

neck_3_L = (1/(1+exp(3.2269-1.9688*maxNIJ_L))); %NIJ calculation

compress_L = abs(min(L_Z_force));
tension_L = max(L_Z_force);

Neck_compress_L = (1/(1+exp(10.9745-2.375*compress_L)));
Neck_tesion_L = (1/(1+exp(10.9745-2.375*tension_L)));


plot(v11660_011(:,1),L_Y_moment)
hold on
plot(v11660_011(:,1),L_Z_force)
plot(v11660_011(:,1),NIJ)
grid on
legend('Moments (N-m)','Forces (N)','NIJ values')
xlabel('Times (seconds)')
title('Neck NIJ, Moments, and Forces of driver of Volkswagen ID.4')
hold off

%% 1b passenger neck
load v11660_055.txt
R_Z_force = v11660_055;
load v11660_057.txt
R_Y_moment = v11660_057;

R_Z_force = j211filtfilt(600,2900,R_Z_force(:,2));
R_Y_moment = j211filtfilt(600,2900,R_Y_moment(:,2));

R_Z_force = R_Z_force/1000;

NIJ_R = zeros(length(R_Y_moment),1);

check_maxNIJ = 0;
NIJ_ind = 0;
load_case_R = '';

for i = 1:length(R_Z_force)
    Fz = R_Z_force(i);
    My = R_Y_moment(i);

    if Fz > 0
        CVForce = 4287/1000; %tension
    else
        CVForce = 3880/1000; %compression
    end
    if My > 0
        CVMoment = 155; %flexion
    else
        CVMoment = 67; %extension
    end

    NIJ_R(i) = abs(Fz/CVForce) + abs(My/CVMoment);

    if check_maxNIJ < NIJ_R(i)
        check_maxNIJ = NIJ_R(i);
        NIJ_ind = i;
        if CVForce == 4.287 && CVMoment == 155
            load_case_R = 'tension/flexion';
        elseif CVForce == 3.880 && CVMoment == 155
            load_case_R = 'compression/flexion';
        elseif CVForce == 4.287 && CVMoment == 67
            load_case_R = 'tension/extension';
        else
            load_case_R = 'compression/extension';
        end
    end
end

maxNIJ_R = max(NIJ_R);

neck_3_R = (1/(1+exp(3.2269-1.9688*maxNIJ_R))); %NIJ calculation

compress_R = abs(min(R_Z_force));
tension_R = max(R_Z_force);

Neck_compress_R = (1/(1+exp(10.958-3.77*compress_R)));
Neck_tesion_R = (1/(1+exp(10.958-3.77*tension_R)));

plot(v11660_011(:,1),R_Y_moment)
hold on
plot(v11660_011(:,1),R_Z_force)
plot(v11660_011(:,1),NIJ_R)
grid on
legend('Moments (N-m)','Forces (N)','NIJ values')
xlabel('Times (seconds)')
title('Neck NIJ, Moments, and Forces of passenger of Volkswagen ID.4')
hold off

%% 1c
load v11660_019.txt %left chest displacement tranducer
load v11660_065.txt %right chest

left_chest = v11660_019;
right_chest = v11660_065;

left_chest = j211filtfilt(600,2900,left_chest(:,2));
right_chest = j211filtfilt(600,2900,right_chest(:,2));

peak_chest_disp_L = abs(min(left_chest));
peak_chest_disp_R = abs(min(right_chest));

L_chest_AIS = 1/(1 + (exp(10.5456-1.568*peak_chest_disp_L^0.4612)));
R_chest_AIS = 1/(1 + (exp(10.5456-1.7212*peak_chest_disp_R^0.4612)));

%% 1d
load v11660_023.txt
Femur_LL = v11660_023; %left passenger left femur
Femur_LL = j211filtfilt(600,2900,Femur_LL(:,2));
Femur_LL = Femur_LL/1000;

load v11660_024.txt
Femur_LR = v11660_024; %left passenger right femur
Femur_LR = j211filtfilt(600,2900,Femur_LR(:,2));
Femur_LR = Femur_LR/1000;

load v11660_069.txt
Femur_RL = v11660_069; %right passenger left femur
Femur_RL = j211filtfilt(600,2900,Femur_RL(:,2));
Femur_RL = Femur_RL/1000;

load v11660_070.txt
Femur_RR = v11660_070; %right passenger right femur
Femur_RR = j211filtfilt(600,2900,Femur_RR(:,2));
Femur_RR = Femur_RR/1000;

max_compressive_LL = max(abs(Femur_LL));
max_compressive_LR = max(abs(Femur_LR));

max_compressive_RL = max(abs(Femur_RL));
max_compressive_RR = max(abs(Femur_RR));

LL_ind = find(max_compressive_LL==abs(Femur_LL));
LR_ind = find(max_compressive_LR==abs(Femur_LR));

RL_ind = find(max_compressive_RL==abs(Femur_RL));
RR_ind = find(max_compressive_RR==abs(Femur_RR));

LL_FemurAIS = 1/(1 + (exp(5.795-0.5196*max_compressive_LL)));
LR_FemurAIS = 1/(1 + (exp(5.795-0.5196*max_compressive_LR)));

RL_FemurAIS = 1/(1 + (exp(5.7949-0.7619*max_compressive_RL)));
RR_FemurAIS = 1/(1 + (exp(5.7949-0.7619*max_compressive_RR)));

plot(v11660_070(:,1),Femur_LL)
hold on 
plot(v11660_070(:,1),Femur_LR)
plot(v11660_070(:,1),Femur_RL)
plot(v11660_070(:,1),Femur_RR)
grid on
xlabel('Times (seconds)')
ylabel('Force (kN)')
title('Forces onto the Femur in ID.4 Driver')
legend('Driver Left Femur','Driver Right Femur','Passenger Left Femur','Passenger Right Femur')
title('Forces onto both Femurs in Honda Occupants')
hold off

plot(v11660_070(:,1),Femur_RL)
hold on 
plot(v11660_070(:,1),Femur_RR)
grid on
xlabel('Times (seconds)')
ylabel('Force (kN)')
title('Forces onto the Femur in ID.4 Passenger')
legend('Left Femur','Right Femur')
hold off


%%
maxfemur_L = max(LL_FemurAIS,LR_FemurAIS);
maxNeck_L = max([Neck_compress_L, Neck_tesion_L, neck_3_L]);

Pinjury_L = 1 - (1 - VW_L_head)*(1 - L_chest_AIS)*(1 - maxfemur_L)*(1 - maxNeck_L);
RR_L = Pinjury_L/0.15;

maxfemur_R = max(RL_FemurAIS,RR_FemurAIS);
maxNeck_R = max([Neck_compress_R, Neck_tesion_R, neck_3_R]);

Pinjury_R = 1 - (1 - VW_R_head)*(1 - R_chest_AIS)*(1 - maxfemur_R)*(1 - maxNeck_R);
RR_R = Pinjury_R/0.15;

%% 2a
load v11660_013.txt
chest_L = v11660_013;

chest_L = j211filtfilt(600,4000,chest_L(:,2));

chest_clips_L = [];

count = 0;
for i = 1:3:length(chest_L)-2
    count = count + 1;
    clip = [chest_L(i) chest_L(i + 1) chest_L(i + 2)];
    chest_clips_L = vertcat(chest_clips_L,max(clip));
end


chest_accel_L = abs(min(chest_clips_L));
CTI = chest_accel_L/90 + peak_chest_disp_L/103;

AIS_Ac = 1/(1+exp(3.1493-0.063*chest_accel_L));
AIS_CTI = 1/(1+exp(8.224-7.125*CTI));

%% 2a
Age = [20 35 50 70 90];

for i = 1:length(Age)
    Pthorax(i) = 1/(1 + exp(-(-12.597 + 0.05861*Age(i) + 1.568*(peak_chest_disp_L^0.4612))));
end


