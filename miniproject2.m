%% Total Barrier Force of Honda in Newtons
Force = zeros(3312,1); %preallocate 0 vector for force
TBF = zeros(3312,2); %preallocate 0 matrix for total barrier force of VW
HMass = 1438;
v0_kph = 56.4;
vf_kph = 7.66256959;

% loop through needed txt files to add forces to total barrier force
for x=97:1:132
    if x <= 99
        filename=strjoin(["v05573_0",x,".txt"],"");
    else
        filename=strjoin(["v05573_",x,".txt"],""); %make file names
    end
    temp_H = load (filename);
    temp_H = j211filtfilt(60,10000,temp_H); %figure out what the deal with the filter is 60 or 180
    Force = Force + temp_H(:,2); %take force column of each file and add to Force vector
end

TBF_H = cat(2,[temp_H(:,1),Force]); %concat time column with force column
TBF_H = TBF_H(304:2951,:); %remove first 304 rows which is before car contacts barrier

IntegrateForce = cumtrapz(TBF_H(:,1),TBF_H(:,2)); %find Impulse

maxImpulse = min(IntegrateForce);
correctionFactor_H = (HMass*(vf_kph-v0_kph)/3.6)/maxImpulse; %correction factor
TBF_H(:,2) = TBF_H(:,2)*correctionFactor_H;

plot(TBF_H(:,1),TBF_H(:,2)) %time vs Force
xlabel('Time (s)')
ylabel('Force (N)')
grid on

min(TBF_H(:,2))

%% Impulse for VW in N-s
plot(TBF_H(:,1),IntegrateForce) %time vs impulse
xlabel('Time (s)')
ylabel('Impulse (N-s)')
grid on

timeVSimpulse = [TBF_H(:,1) IntegrateForce];
maxForce = min(TBF_H(:,2));

%% Obtain Acceleration of Occupant Comp and Integrate acceleration to get velocity
load v05573_089.txt %occupant comp
OccComp=v05573_089(300:2947,:); %remove first 500 rows before crash
%OccComp(:,2) = OccComp(:,2);%*9.8 % convert Gs to m/s^2
OccCompFilt = j211filtfilt(180,10000,OccComp); %filter data
v0 = 35; %mph to add to velocity values

plot(OccComp(:,1),OccCompFilt(:,2)) %time vs accel
xlabel('Time (s)')
ylabel('Acceleration (g)')
grid on

OccCompInt = OccCompFilt;
OccCompInt(:,2) = OccCompFilt(:,2)*9.8;

IntegrateAccel = cumtrapz(OccCompInt(:,1),OccCompInt(:,2)); %find accel in m/s
IntegrateAccel = IntegrateAccel*2.237 + v0; %convert to mph and add initial velo

plot(OccComp(:,1),IntegrateAccel) % velo vs time
xlabel('Time (s)')
ylabel('Velocity (mph)')
grid on

veloVstime_H = [OccComp(:,1) IntegrateAccel];


%% Integrate velocity to get displacement
GetReadyToIntegrate = IntegrateAccel/3600; % convert to miles per second before integrating

IntegrateVelo = cumtrapz(OccComp(:,1),GetReadyToIntegrate); %Find displacement
IntegrateVelo = IntegrateVelo*63360; %convert miles to inches

plot(OccComp(:,1),IntegrateVelo) %time vs Displacement
xlabel('Time (s)')
ylabel('Displacement (in)')
grid on

timeVsDisplacement_H = [OccComp(:,1) IntegrateVelo];

plot(IntegrateVelo,IntegrateAccel) % Displacement vs Velo
xlabel('Displacement (in)')
ylabel('Velocity (mph)')
grid on

dispVsVelo_H = [ IntegrateVelo IntegrateAccel ]; % for combined graph at end

plot(IntegrateVelo,OccCompFilt) % Displacement vs Acceleration
xlabel('Displacement (in)')
ylabel('Acceleration (g)')
grid on

H_force_from_accel = OccCompInt(:,2)*HMass; % get Force by F=MA

plot(IntegrateVelo,(H_force_from_accel*-1)/1000) % Displacement vs Force
xlabel('Displacement (in)')
ylabel('Force (kN)')
grid on

plot(IntegrateVelo,(TBF_H(:,2)*-1)/1000) %Displacement vs Barrier Force
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on

TBF_H_lbs = TBF_H(:,2)*0.224809; %convert N to pounds of force
forceVsDisplacement_H = [IntegrateVelo TBF_H_lbs]; %don't use for now

plot(IntegrateVelo, TBF_H_lbs)
xlabel('Displacement (in)')
ylabel('Force (lbs)')
grid on


MaxDynCrush = max(IntegrateVelo); %max displacement 
dispVSforce = [IntegrateVelo TBF_H(:,2)*-1/1000];

min(TBF_H(:,2))

%----------------------------------------------------------
%% Total Barrier Force of VW in Newtons
Force = zeros(3500,1); %preallocate 0 vector for force
TBF_VW = zeros(3500,2); %preallocate 0 matrix for total barrier force of VW
VWMass = 2305;
v0_kph = 56.4;
vf_kph = -13.2;

% loop through needed txt files to add forces to total barrier force
for x=107:3:634
    filename=strjoin(["v11660_",x,".txt"],""); %make file names
    temp = load (filename);
    temp = j211filtfilt(60,10000,temp); %figure out what the deal with the filter is 60 or 180
    Force = Force + temp(:,2); %take force column of each file and add to Force vector
end

TBF_VW = cat(2,[temp(:,1),Force]); %concat time column with force column
TBF_VW = TBF_VW(500:3500,:); %remove first 500 rows which is before car contacts barrier

%TBF_VW_Filt = j211filtfilt(60,10000,TBF_VW); %filter
IntegrateForce = cumtrapz(TBF_VW(:,1),TBF_VW(:,2));

maxImpulse = min(IntegrateForce);
correctionFactor_VW = (VWMass*(vf_kph-v0_kph)/3.6)/maxImpulse;

TBF_VW(:,2) = TBF_VW(:,2)*correctionFactor_VW;
maxForce = min(TBF_VW(:,2));

plot(TBF_VW(:,1),TBF_VW(:,2)) %time vs Force
xlabel('Time (s)')
ylabel('Force (N)')
grid on


%% Impulse for VW in N-s
 %find Impulse
plot(TBF_VW(:,1),IntegrateForce) %time vs impulse
xlabel('Time (s)')
ylabel('Impulse (N-s)')
grid on

timeVSimpulse = [TBF_VW(:,1) IntegrateForce];

%% Obtain Acceleration of Occupant Comp and Integrate acceleration to get velocity
load v11660_093.txt %occupant comp
%load v11660_099.txt
OccComp=v11660_093(500:3500,:); %remove first 500 rows before crash
%OccComp(:,2) = OccComp(:,2);%*9.8 % convert Gs to m/s^2
OccCompFilt = j211filtfilt(180,10000,OccComp); %filter data
v0 = 35; %m/s to add to velocity values

plot(OccComp(:,1),OccCompFilt(:,2)) %time vs accel
xlabel('Time (s)')
ylabel('Acceleration (g)')
grid on

OccCompInt = OccCompFilt;
OccCompInt(:,2) = OccCompFilt(:,2)*9.8;

IntegrateAccel = cumtrapz(OccCompInt(:,1),OccCompInt(:,2)); %find accel in m/s
IntegrateAccel = IntegrateAccel*2.237 + v0; %convert to mph and add initial velo

plot(OccComp(:,1),IntegrateAccel) % velo vs time
xlabel('Time (s)')
ylabel('Velocity (mph)')
grid on

veloVstime_VW = [OccComp(:,1) IntegrateAccel];


%% Integrate velocity to get displacement
GetReadyToIntegrate = IntegrateAccel/3600; % convert to miles per second before integrating

IntegrateVelo = cumtrapz(OccComp(:,1),GetReadyToIntegrate); %Find displacement
IntegrateVelo = IntegrateVelo*63360; %convert miles to inches

plot(OccComp(:,1),IntegrateVelo) %time vs Displacement
xlabel('Time (s)')
ylabel('Displacement (in)')
grid on

timeVsDisplacement_VW = [OccComp(:,1) IntegrateVelo];

plot(IntegrateVelo,IntegrateAccel) % Displacement vs Velo
xlabel('Displacement (in)')
ylabel('Velocity (mph)')
grid on

dispVsVelo_VW = [ IntegrateVelo IntegrateAccel ]; % for combined graph at end

plot(IntegrateVelo,OccCompFilt) % Displacement vs Acceleration
xlabel('Displacement (in)')
ylabel('Acceleration (g)')
grid on

VW_force_from_accel = OccCompInt(:,2)*VWMass; % get Force by F=MA

plot(IntegrateVelo,(VW_force_from_accel*-1)/1000) % Displacement vs Force
xlabel('Displacement (in)')
ylabel('Force (kN)')
grid on

plot(IntegrateVelo,(TBF_VW(:,2)*-1)/1000) %Displacement vs Barrier Force
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on

TBF_VW_lbs = TBF_VW(:,2)*0.224809; %convert N to pounds of force
forceVsDisplacement_VW = [IntegrateVelo TBF_VW_lbs];

plot(IntegrateVelo, TBF_VW_lbs)
xlabel('Displacement (in)')
ylabel('Force (lbs)')
grid on

MaxDynCrush = max(IntegrateVelo); %max displacement 
dispVSforce = [IntegrateVelo (VW_force_from_accel*-1)/1000];



%% comparing honda and VW on same plots

plot(forceVsDisplacement_VW(:,1),forceVsDisplacement_VW(:,2),'r')
xlabel('Displacement (in)')
ylabel('Force (lbs)')
grid on
title('Honda Vs VW Displacement and Fore')

hold on

plot(forceVsDisplacement_H(:,1),forceVsDisplacement_H(:,2),'k')
legend('VW','Honda')


plot(dispVsVelo_H(:,1),dispVsVelo_H(:,2),'k') %comparing honda and VW for disp vs velo
xlabel('Displacement (in)')
ylabel('Velocity (mph)')
title('Displacement versus Velocity of Honda and Volkswagen')
grid on

hold on

plot(dispVsVelo_VW(:,1),dispVsVelo_VW(:,2),'r')
legend('Honda','VW')

% compare H vs VW time vs displacement
plot(timeVsDisplacement_H(:,1),timeVsDisplacement_H(:,2),'k')
xlabel('Time (seconds)')
ylabel('Displacement (in)')
title('Time versus Displacement')
grid on
hold on
plot(timeVsDisplacement_VW(:,1),timeVsDisplacement_VW(:,2),'r')
xlim([0 .3])
ylim([0 35])
legend('Honda','VW')

%compare H vs VW time vs velocity
plot(veloVstime_H(:,1),veloVstime_H(:,2),'k')
xlabel('Time (seconds)')
ylabel('Velocity (mph)')
title('Time versus Velocity')
grid on
hold on
plot(veloVstime_VW(:,1),veloVstime_VW(:,2),'r')
legend('Honda','VW')

%accel for Honda vs VW
OccCompFilt = OccCompFilt*9.8;
OccCompFilt_VW = OccCompFilt_VW*9.8;


plot(OccComp(:,1),OccCompFilt(:,2)) %time vs accel
hold on 
plot(OccComp_VW(:,1),OccCompFilt_VW(:,2))
title('Acceleration versus Velocity')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
xlim([0 .3])
ylim([-400 10])
legend('Honda','VW')
grid on


