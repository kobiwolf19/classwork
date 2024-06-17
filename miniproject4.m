Ku = 18028306; %unloading stiffness kg/s^2
Kl = 993462; %loading stiffness kg/s^2
m = 2305; % vehicle mass kg
v0 = 1000*56.4/3600; % initial speed m/s

Wu = sqrt((Ku/m)); % 1/s
Wl = sqrt((Kl/m)); % 1/s

accelAtMC = -v0*Wl*sin(Wl*(pi/(2*Wl))); %m/s^2

maxForce = abs(m*accelAtMC); % x xMax (N)

dispMC = (v0/Wl)*sin(Wl*(pi/(2*Wl)))*39.3701; % y Fmax (inches)

x0 = ((v0/Wl) - ((Kl/Ku)*(v0/Wl)))*39.3701; % x0 position (inches)


%% Honda
m_H = 1438; %kg
Ku_H = 36944855; %unloading stiffness kg/s^2
Kl_H = 692839; %loading stiffness kg/s^2

Wu_H = sqrt((Ku_H/m_H)); % 1/s
Wl_H = sqrt((Kl_H/m_H)); % 1/s

accelAtMC_H = -v0*Wl_H*sin(Wl_H*(pi/(2*Wl_H))); %m/s^2

maxForce_H = abs(m_H*accelAtMC_H); % x xMax (N)

dispMC_H = (v0/Wl_H)*sin(Wl_H*(pi/(2*Wl_H)))*39.3701; % y Fmax (inches)

x0_H = ((v0/Wl_H) - ((Kl_H/Ku_H)*(v0/Wl_H)))*39.3701; % x0 position (inches)


plot([0 dispMC], [0 maxForce])
hold on
plot([dispMC x0], [maxForce 0])
plot([0 dispMC_H], [0 maxForce_H])
plot([dispMC_H x0_H], [maxForce_H 0])
xlabel('Displacment (inches)')
ylabel('Force (N)')
title('F-D of Loading and Unloading Phases')
legend('VW Loading','VW Unloading', 'Honda Loading','Honda Unloading')
grid on
hold off

%% part b

% displacments
v0 = (40/3600)*63360; %inches per second
Vr = -(Kl/Ku)*(v0/Wl)*Wu*sin(Wu*pi/(2*Wu));
Vr_H = -(Kl_H/Ku_H)*(v0/Wl_H)*Wu_H*sin(Wu_H*pi/(2*Wu_H));

fplot(@(t) (v0/Wl)*sin(Wl*t), [0 pi/(2*Wl)])
hold on
fplot(@(t) (Kl/Ku)*(v0/Wl)*cos(Wu*(t-(pi/(2*Wl)))) + x0, [pi/(2*Wl) (pi/(2*Wl)+pi/(2*Wu))])
fplot(@(t) (Vr*(t-(pi/(2*Wl)+pi/(2*Wu))) + x0), [(pi/(2*Wl)+pi/(2*Wu)) 0.3])
fplot(@(t) (v0/Wl_H)*sin(Wl_H*t), [0 pi/(2*Wl_H)])
fplot(@(t) (Kl_H/Ku_H)*(v0/Wl_H)*cos(Wu_H*(t-(pi/(2*Wl_H)))) + x0_H, [pi/(2*Wl_H) (pi/(2*Wl_H)+pi/(2*Wu_H))])
fplot(@(t) (Vr_H*(t-(pi/(2*Wl_H)+pi/(2*Wu_H))) + x0_H), [(pi/(2*Wl_H)+pi/(2*Wu_H)) 0.3])
xlim([0 .3])
ylim([0 35])
grid on
xlabel('Time (s)')
ylabel('Displacement (in)')
title('Displacement during Frontal Impact Test at 40 mph')
legend('VW Loading', 'VW Unloading', 'VW Separation', 'Honda Loading', 'Honda Unloading', 'Honda Separation')
hold off


% velocity 
v0 = 35; %mph
Vr = -(Kl/Ku)*(v0/Wl)*Wu*sin(Wu*pi/(2*Wu));
Vr_H = -(Kl_H/Ku_H)*(v0/Wl_H)*Wu_H*sin(Wu_H*pi/(2*Wu_H));

fplot(@(t) v0*cos(Wl*t), [0 pi/(2*Wl)])
hold on
fplot(@(t) -(Kl/Ku)*(v0/Wl)*Wu*sin(Wu*(t-(pi/(2*Wl)))), [pi/(2*Wl) (pi/(2*Wl)+pi/(2*Wu))])
fplot(@(t) Vr + 0*t, [(pi/(2*Wl)+pi/(2*Wu)) 0.3])
fplot(@(t) v0*cos(Wl_H*t), [0 pi/(2*Wl_H)])
fplot(@(t) -(Kl_H/Ku_H)*(v0/Wl_H)*Wu_H*sin(Wu_H*(t-(pi/(2*Wl_H)))), [pi/(2*Wl_H) (pi/(2*Wl_H)+pi/(2*Wu_H))])
fplot(@(t) Vr_H + 0*t, [(pi/(2*Wl_H)+pi/(2*Wu_H)) 0.3])
xlim([0 .3])
ylim([-10 40])
grid on
xlabel('Time (s)')
ylabel('Velocity (mph)')
title('Velocity during Frontal Impact Test at 40 mph')
legend('VW Loading', 'VW Unloading', 'VW Separation', 'Honda Loading', 'Honda Unloading', 'Honda Separation')
hold off


% acceleration
v0 = 35/2.237; %m/s

fplot(@(t) -v0*Wl*sin(Wl*t), [0 pi/(2*Wl)])
hold on
fplot(@(t) -(Kl/Ku)*(v0/Wl)*Wu^2*cos(Wu*(t-(pi/(2*Wl)))), [pi/(2*Wl) (pi/(2*Wl)+pi/(2*Wu))])
fplot(@(t) 0*t, [(pi/(2*Wl)+pi/(2*Wu)) 0.3])
fplot(@(t) -v0*Wl_H*sin(Wl_H*t), [0 pi/(2*Wl_H)])
fplot(@(t) -(Kl_H/Ku_H)*(v0/Wl_H)*Wu_H^2*cos(Wu_H*(t-(pi/(2*Wl_H)))), [pi/(2*Wl_H) (pi/(2*Wl_H)+pi/(2*Wu_H))])
fplot(@(t) 0*t, [(pi/(2*Wl_H)+pi/(2*Wu_H)) 0.3])
xlim([0 .3])
ylim([-400 10])
grid on
xlabel('Time (s)')
ylabel('Accleration (m/s^2)')
title('Acceleration during Frontal Impact Test at 40 mph')
legend('VW Loading', 'VW Unloading', 'VW Separation', 'Honda Loading', 'Honda Unloading', 'Honda Separation')
hold off

