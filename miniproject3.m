load v11660_093.txt %occupant comp

OccComp=v11660_093(500:3500,:); %after crash
OccCompFilt = j211filtfilt(180,10000,OccComp); %filter data

v0 = 35; %mph to add to velocity values
v0_kph = 56.4; %initial kph
m_VW = 2305;

OccCompInt = [OccComp(:,1) OccCompFilt(:,2)*9.8];

IntegrateAccel_OC = cumtrapz(OccCompInt(:,1),OccCompInt(:,2)); %obtain velo
IntegrateAccel_OC = IntegrateAccel_OC*2.237 + v0; %convert to mph and add initial velo

plot(OccComp(:,1),IntegrateAccel_OC) % velo vs time of Occ Comp
xlabel('Time (s)')
ylabel('Velocity (mph)')
grid on

%% Integrate velocity to get displacement
GetReadyToIntegrate_OC = IntegrateAccel_OC/3600; % convert to miles per second before integrating

IntegrateVelo_OC = cumtrapz(OccComp(:,1),GetReadyToIntegrate_OC); %Find displacement
IntegrateVelo_OC = IntegrateVelo_OC*63360; %convert miles to inches

%% chest calculations
load v11660_013.txt %chest accel

ChestAccel = v11660_013(500:3500,:);
ChestAccelFilt = j211filtfilt(180,10000,ChestAccel);

ChestAccelInt = [ChestAccel(:,1) ChestAccelFilt(:,2)*9.8];

IntegrateChestAccel = cumtrapz(ChestAccelInt(:,1),ChestAccelInt(:,2)); %chest Velo for inches
IntegrateChestAccel = IntegrateChestAccel*2.237 + v0; %convert chest velo to mph with initial velo

GetReadyToIntegrate_Chest = IntegrateChestAccel/3600; % convert to miles per second before integrating

IntegrateChestVelo = cumtrapz(ChestAccel(:,1),GetReadyToIntegrate_Chest); % chest displacement calc
IntegrateChestVelo = IntegrateChestVelo*63360; %convert miles to inches

IntegrateAccelChest_GsOC = cumtrapz(IntegrateVelo_OC,ChestAccelFilt(:,2)); %integrate accel in terms of G force and inches

%% calculations/plots
energyInitial = 1/2*(v0^2);

plot(IntegrateVelo_OC,IntegrateAccelChest_GsOC); % 2.a.2 vehicle energy absorbed vs disaplacement of vehicle
xlabel('Displacement of Vehicle (inches)')
ylabel('Energy Absorbed (G-in)')
grid on
title('Energy Absorbed by vehicle for VW ID.4')

vehicleChest =  IntegrateChestVelo - IntegrateVelo_OC; %displacement of occupant WRT vehicle in inches

integrateRestraints = cumtrapz(vehicleChest,ChestAccelFilt(:,2)); %integrate chest accel vs restratint/occupant displacement

plot(vehicleChest,integrateRestraints); % 2.a.1
xlabel('Displacement of Chest/Restrains (inches)')
ylabel('Energy Absorbed (G-in)')
grid on
title('Energy Absorbed by restrains/occupant for the driver of the VW ID.4')

%% restraint quotient
plot(OccComp(:,1),IntegrateAccel_OC)
xlabel('Velocity (mph)')
ylabel('Time (s)')
grid on
hold on
plot(OccComp(:,1),IntegrateChestAccel)
legend('Occupant Compartment','Driver Chest')

VC = abs(IntegrateAccel_OC - IntegrateChestAccel); %difference of velocities of chest and vehicle
RQ = VC/v0; %Restraint quotient equation

plot(OccComp(:,1),RQ); %2c RQ vs time
xlabel('Time (s)')
ylabel('RQ')
grid on

RKEF = (VC.^2)/25;

plot(OccComp(:,1),RKEF); %2c RKEF vs time
xlabel('Time (s)')
ylabel('RKEF')
grid on

maxRKEF = max(RKEF);
maxRQ = max(RQ);

%% right side passanger
load v11660_059.txt

ChestAccel_R = v11660_059(500:3500,:);
ChestAccelFilt_R = j211filtfilt(180,10000,ChestAccel_R);

ChestAccelInt_R = [ChestAccel_R(:,1) ChestAccelFilt_R(:,2)*9.8];

IntegrateChestAccel_R = cumtrapz(ChestAccelInt_R(:,1),ChestAccelInt_R(:,2)); %chest Velo for inches
IntegrateChestAccel_R = IntegrateChestAccel_R*2.237 + v0; %convert chest velo to mph with initial velo

GetReadyToIntegrate_Chest_R = IntegrateChestAccel_R/3600; % convert to miles per second before integrating

IntegrateChestVelo_R = cumtrapz(ChestAccel_R(:,1),GetReadyToIntegrate_Chest_R); % chest displacement calc
IntegrateChestVelo_R = IntegrateChestVelo_R*63360; %convert miles to inches

IntegrateAccelChest_GsOC_R = cumtrapz(IntegrateVelo_OC,ChestAccelFilt_R(:,2)); %integrate accel in terms of G force and inches

%% Calculations/plots for right passenger

plot(IntegrateVelo_OC,IntegrateAccelChest_GsOC_R); % 2.a.2 vehicle energy absorbed vs disaplacement of vehicle (Right passenger)
xlabel('Displacement of Vehicle (inches)')
ylabel('Energy Absorbed (G-in)')
grid on
title('Energy Absorbed by vehicle')

vehicleChest_R =  IntegrateChestVelo_R - IntegrateVelo_OC; %displacement of occupant WRT vehicle in inches

integrateRestraints_R = cumtrapz(vehicleChest_R,ChestAccelFilt_R(:,2)); %integrate chest accel vs restratint/occupant displacement

plot(vehicleChest_R,integrateRestraints_R); % 2.a.1 right passenger
xlabel('Displacement of Chest/Restrains (inches)')
ylabel('Energy Absorbed (G-in)')
grid on
title('Energy Absorbed by restrains/occupant')

%% RQ and RKEF right passenger
VC_R = abs(IntegrateAccel_OC - IntegrateChestAccel_R); %difference of velocities of chest and vehicle
RQ_R = VC_R/v0; %Restraint quotient equation

plot(OccComp(:,1),RQ_R); %2c RQ vs time
xlabel('Time (s)')
ylabel('RQ')
grid on

RKEF_R = (VC_R.^2)/25;

plot(OccComp(:,1),RKEF_R); %2c RKEF vs time
xlabel('Time (s)')
ylabel('RKEF')
grid on

maxRKEF_R = max(RKEF_R);
maxRQ_R = max(RQ_R);

%% unbelted plots/calculations
IntegrateAccel_OC = IntegrateAccel_OC - v0; % to use for wrt values

GetReadyToIntegrate_OC = IntegrateAccel_OC/3600; % convert to miles per second before integrating

IntegrateVelo_OC = cumtrapz(OccComp(:,1),GetReadyToIntegrate_OC); %Find displacement
IntegrateVelo_OC = IntegrateVelo_OC*63360;

plot(OccCompFilt(:,1),-OccCompFilt(:,2)) % plot showing a',v',x'
hold on
plot(OccCompFilt(:,1),-IntegrateAccel_OC)
hold on
plot(OccCompFilt(:,1),-IntegrateVelo_OC)
axis([0 0.15 0 50]);
grid on
legend('a'' (G)', 'v'' (mph)','x'' (in)')
xlabel('Time (s)')

plot(-IntegrateVelo_OC,-IntegrateAccel_OC) % driver displacement vs velocity wrt vehicle
axis([0 60 0 50])
xlabel('Driver Displacement w.r.t Vehicle (in)')
ylabel('Driver Velocity w.r.t Vehicle (mph)')
title('Unrestrained Occupant Kinetics for VW ID.4')
grid on

check_15in = [-IntegrateVelo_OC -IntegrateAccel_OC];


%% L vs R graphs and Honda vs VW
plot(vehicleChest,integrateRestraints); % 2.a.1 L vs R
xlabel('Displacement of Chest/Restrains (inches) WRT vehicle')
ylabel('Energy Absorbed (G-in)')
grid on
hold on
plot(vehicleChest_R,integrateRestraints_R)
hold on
plot(vehicleChest_H,integrateRestraints_H);
hold on
plot(vehicleChest_R_H,integrateRestraints_R_H)
legend('VW Driver','VW Right front passenger','Honda Driver','Honda Right front passenger')
title('Energy Absorbed by restrains/occupant')


plot(IntegrateVelo_OC,IntegrateAccelChest_GsOC); % 2.a.2 vehicle energy absorbed vs disaplacement of vehicle
xlabel('Displacement of Vehicle (inches)')
ylabel('Energy Absorbed (G-in)')
grid on
hold on
plot(IntegrateVelo_OC,IntegrateAccelChest_GsOC_R);
hold on
plot(IntegrateVelo_OC_H,IntegrateAccelChest_GsOC_H);
hold on
plot(IntegrateVelo_OC_H,IntegrateAccelChest_GsOC_R_H);
legend('VW Driver side','VW Right front passenger side','Honda Driver side','Honda Right front passenger side')
title('Energy Absorbed by vehicle structure')


plot(OccComp(:,1),RQ); %2c RQ vs time
xlabel('Time (s)')
ylabel('RQ')
hold on
plot(OccComp(:,1),RQ_R);
hold on
plot(OccComp_H(:,1),RQ_H);
hold on
plot(OccComp_H(:,1),RQ_R_H);
legend('VW Driver','VW Right Front Passenger','Honda Driver','Honda Right Front Passenger')
title('Restraint Quotient')
grid on


plot(OccComp(:,1),RKEF); %2c RKEF vs time
xlabel('Time (s)')
ylabel('RKEF')
hold on
plot(OccComp(:,1),RKEF_R);
hold on
plot(OccComp_H(:,1),RKEF_H);
hold on 
plot(OccComp_H(:,1),RKEF_R_H);
legend('VW Driver','VW Right Front Passenger','Honda Driver','Honda Right Front Passenger')
title('Relative Kinetic Energy Factor')
grid on




