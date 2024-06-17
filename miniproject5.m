%35 mph
global xmax
xmax = 0;
x0 = 0;
v0 = 35;
v0 = v0/3600*63360; %in/s
[t35,y35] = ode23(@onespring, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t35 y35(:,2)/17.6 y35(:,1)];

for i = 2:length(y35)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration35 = COV./COT;
Acceleration35 = Acceleration35*0.045585;

plot(t35(2:538),Acceleration35) %plot acceleration in Gs


plot(t35,(y35(:,2))/17.6)
plot(t35,y35(:,1))

%% 30 mph
global xmax
xmax = 0;
x0 = 0;
v0 = 30;
v0 = v0/3600*63360; %in/s
[t30,y30] = ode23(@onespring, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t30 y30(:,2)/17.6 y30(:,1)];

COV = zeros(1,length(y30));
COT = zeros(1,length(y30));

for i = 2:length(y30)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration30 = COV./COT;
Acceleration30 = Acceleration30*0.045585;

plot(t30(1:529),Acceleration30)

%% 25 mph
global xmax
xmax = 0;
x0 = 0;
v0 = 25;
v0 = v0/3600*63360; %in/s
[t25,y25] = ode23(@onespring, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t25 y25(:,2)/17.6 y25(:,1)];

COV = zeros(1,length(y25));
COT = zeros(1,length(y25));

for i = 2:length(y25)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration25 = COV./COT;
Acceleration25 = Acceleration25*0.045585;

plot(t25(1:516),Acceleration25)

%all 3 VW accelerations
plot(t35(2:538),Acceleration35)
hold on
plot(t30(1:529),Acceleration30)
plot(t25(1:516),Acceleration25)
hold off
grid on 
ylabel('Acceleration (G''s)')
xlabel('Time (seconds)')
title('Calculated Accelerations of VW ID.4 using ODE23')
legend('35 mph','30 mph','25 mph')
xlim([0 0.1])

%VW velocities
plot(t35,(y35(:,2))/17.6)
hold on
plot(t30,(y30(:,2))/17.6)
plot(t25,(y25(:,2))/17.6)
hold off
grid on 
ylabel('Velocity (mph)')
xlabel('Time (seconds)')
title('Calculated Velocities of VW ID.4 using ODE23')
legend('35 mph','30 mph','25 mph')
xlim([0 0.2])

%VW displacements
plot(t35,y35(:,1))
hold on
plot(t30,y30(:,1))
plot(t25,y25(:,1))
hold off
grid on 
ylabel('Displacement (in)')
xlabel('Time (seconds)')
title('Calculated Displacements of VW ID.4 using ODE23')
legend('35 mph','30 mph','25 mph')

v0 = 25;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position25(p,1) = i;
    position25(p,2) = v0*i;
end

plot(position25(:,1),position25(:,2))
hold on
plot(t25,y25(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 25 mph in VW ID.4')
xlim([0 0.15])



% 30 passenger
v0 = 30;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position30(p,1) = i;
    position30(p,2) = v0*i;
end

plot(position30(:,1),position30(:,2))
hold on
plot(t30,y30(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 30 mph in VW ID.4')
xlim([0 0.15])

% 35 passenger
v0 = 35;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position35(p,1) = i;
    position35(p,2) = v0*i;
end

plot(position35(:,1),position35(:,2))
hold on
plot(t35,y35(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 35 mph in VW ID.4')
xlim([0 0.15])


%% DRLR sensor correction to find a,v,x

load v11660_093.txt %occupant comp

OccComp=v11660_093(500:3500,:); %remove first 500 rows before crash

OccCompFilt = j211filtfilt(180,10000,OccComp); %filter data

OccCompFilt(:,2) = OccCompFilt(:,2) * 1.053;
v0 = 35;

OccCompInt = OccCompFilt;
OccCompInt(:,2) = OccCompFilt(:,2)*9.8;

VW_Velo = cumtrapz(OccCompInt(:,1),OccCompInt(:,2)); %find velo in m/s
VW_Velo = VW_Velo*2.237 + v0;

GetReadyToIntegrate = VW_Velo/3600; % convert to miles per second before integrating

VW_Disp = cumtrapz(OccComp(:,1),GetReadyToIntegrate); %Find displacement
VW_Disp = VW_Disp*63360;

plot(OccComp(:,1),VW_Disp)

%% Barrier force CF
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

plot(VW_Disp,(TBF_VW(:,2)*-1)/1000) %Displacement vs Barrier Force
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on
ylabel('Force (kN)')
xlabel('Displacement (in)')
title('VW ID.4 Force Deflection Curve')



%% 

VW_FD_load = [0 0; 2.16 22.83; 6.34 357.01; 7.56 348.11; 12.11 524.01; 15.57 445.1; 19.04 532.42;
    19.85 553.76; 20.77 601.6; 21.92 641.34; 22.91 667.38; 25.36 346.43; 27.78 449.58; 28.53 587.73];

VW_impulse = trapz(VW_FD_load(:,1),VW_FD_load(:,2));

VW_impulse = VW_impulse*25.4;

VW_mv2 = 0.5*2305*(15.55556^2);

CFC = VW_impulse/VW_mv2;

VW_FD = [0 0; 2.16 22.83; 6.34 357.01; 7.56 348.11; 12.11 524.01; 15.57 445.1; 19.04 532.42;
    19.85 553.76; 20.77 601.6; 21.92 641.34; 22.91 667.38; 25.36 346.43; 27.78 449.58; 28.53 587.73;
    28.62 489.33; 26.58 124.45; 18.37 3.79; 14.12 -1.19];

VW_FD(:,2) = VW_FD(:,2)*CFC;

plot(VW_FD(:,1),VW_FD(:,2))
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on
title('Approximate VW ID.4 Force Deflection Curve')

%% Honda 35 mph

global xmax
xmax = 0;
x0 = 0;
v0 = 35;
v0 = v0/3600*63360; %in/s
[t35,y35] = ode23(@onespring_H, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t35 y35(:,2)/17.6 y35(:,1)];

for i = 2:length(y35)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration35 = COV./COT;
Acceleration35 = Acceleration35*0.045585;

plot(t35(2:313),Acceleration35*0.045585) %plot acceleration in Gs


plot(t35,(y35(:,2))/17.6)
plot(t35,y35(:,1))

%% 30 mph
global xmax
xmax = 0;
x0 = 0;
v0 = 30;
v0 = v0/3600*63360; %in/s
[t30,y30] = ode23(@onespring_H, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t30 y30(:,2)/17.6 y30(:,1)];

COV = zeros(1,length(y30));
COT = zeros(1,length(y30));

for i = 2:length(y30)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration30 = COV./COT;
Acceleration30 = Acceleration30*0.045585;

plot(t30(1:316),Acceleration30*0.045585)

%% 25 mph
global xmax
xmax = 0;
x0 = 0;
v0 = 25;
v0 = v0/3600*63360; %in/s
[t25,y25] = ode23(@onespring_H, [0 0.300], [x0; v0], odeset('RelTol',1e-8));

tester = [t25 y25(:,2)/17.6 y25(:,1)];

COV = zeros(1,length(y25));
COT = zeros(1,length(y25));

for i = 2:length(y25)
    COV(i-1)=tester(i,2)-tester(i-1,2);
    COT(i-1)=tester(i,1)-tester(i-1,1);
end

Acceleration25 = COV./COT;
Acceleration25 = Acceleration25*0.045585;

plot(t25(1:307),Acceleration25*0.045585)

plot(t35(2:313),Acceleration35*0.045585)
hold on
plot(t30(1:316),Acceleration30*0.045585)
plot(t25(1:307),Acceleration25*0.045585)
hold off
grid on 
ylabel('Acceleration (G''s)')
xlabel('Time (seconds)')
title('Calculated Accelerations of Honda Civic using ODE23')
legend('35 mph','30 mph','25 mph')
xlim([0 0.1])

plot(t35,(y35(:,2))/17.6)
hold on
plot(t30,(y30(:,2))/17.6)
plot(t25,(y25(:,2))/17.6)
hold off
grid on 
ylabel('Velocity (mph)')
xlabel('Time (seconds)')
title('Calculated Velocities of Honda Civic using ODE23')
legend('35 mph','30 mph','25 mph')
xlim([0 0.2])


plot(t35,y35(:,1))
hold on
plot(t30,y30(:,1))
plot(t25,y25(:,1))
hold off
grid on 
ylabel('Displacement (in)')
xlabel('Time (seconds)')
title('Calculated Displacements of Honda Civic using ODE23')
legend('35 mph','30 mph','25 mph')
xlim([0 0.2])


v0 = 25;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position25(p,1) = i;
    position25(p,2) = v0*i;
end

plot(position25(:,1),position25(:,2))
hold on
plot(t25,y25(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 25 mph in Honda Civic')
xlim([0 0.15])

% 30 passenger
v0 = 30;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position30(p,1) = i;
    position30(p,2) = v0*i;
end

plot(position30(:,1),position30(:,2))
hold on
plot(t30,y30(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 30 mph in Honda Civic')
xlim([0 0.15])

% 35 passenger
v0 = 35;
v0 = v0/3600*63360; %in/s

p = 0;
for i = .000:.001:.3
    p = p + 1;
    position35(p,1) = i;
    position35(p,2) = v0*i;
end

plot(position35(:,1),position35(:,2))
hold on
plot(t35,y35(:,1) + 15)
grid on 
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Unbelted Occupant Kinematics at 35 mph in Honda Civic')
xlim([0 0.15])


%% DRLR sensor correction to find a,v,x

load v05573_089.txt %occupant comp
OccComp=v05573_089(300:2947,:); %remove first 500 rows before crash
%OccComp(:,2) = OccComp(:,2);%*9.8 % convert Gs to m/s^2
OccCompFilt = j211filtfilt(180,10000,OccComp); %filter data

OccCompFilt(:,2) = OccCompFilt(:,2) * 0.8651;
v0 = 35; %mph to add to velocity values

OccCompInt = OccCompFilt;
OccCompInt(:,2) = OccCompFilt(:,2)*9.8;

IntegrateAccel = cumtrapz(OccCompInt(:,1),OccCompInt(:,2)); %find accel in m/s
IntegrateAccel = IntegrateAccel*2.237 + v0; %convert to mph and add initial velo

veloVstime_H = [OccComp(:,1) IntegrateAccel];

GetReadyToIntegrate = IntegrateAccel/3600; % convert to miles per second before integrating

IntegrateVelo = cumtrapz(OccComp(:,1),GetReadyToIntegrate); %Find displacement
IntegrateVelo = IntegrateVelo*63360; %convert miles to inches

plot(OccComp(:,1),IntegrateVelo) %time vs Displacement
xlabel('Time (s)')
ylabel('Displacement (in)')
grid on

%% Barrier force CF
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

plot(IntegrateVelo,(TBF_H(:,2)*-1)/1000) %Displacement vs Barrier Force
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on
title('Honda Civic Force Deflection Curve')



%% 

H_FD_load = [0 0; 7.2 171.59; 9.05 166.67; 15.47 417.87; 18.53 262.21; 21.38 214.42; 23.44 174.3;
    25.57 232.52; 27.38 275.87; 28.89 218.47; 30.47 239.52; 31.02 0.43];

H_impulse = trapz(H_FD_load(:,1),H_FD_load(:,2));

H_impulse = H_impulse*25.4;

H_mv2 = 0.5*2305*(15.55556^2);

CFC = H_impulse/H_mv2;


H_FD_load(:,2) = H_FD_load(:,2)*CFC;

plot(H_FD_load(:,1),H_FD_load(:,2))
xlabel('Displacement (in)')
ylabel('Barrier Force (kN)')
grid on
title('Approximate Honda Civic Force Deflection Curve')

stiffness = zeros(12,1);
for i = 2:12 % N/m
    stiffness(i) = H_FD_load(i,2)/H_FD_load(i,1)*39370; %force/displacement
end
