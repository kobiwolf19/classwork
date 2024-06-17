function filtered = j211filter(x,fc,fs)
n = 4; %number of poles
[b,a] = butter(n,fc/(fs/2));%returns transfer function coefficients
filtered = filtfilt(b,a,x); %returns filtered signal
end