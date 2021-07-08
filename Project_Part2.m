%% EE-430 Project Part 2, Berkay Yaldiz - 2232940, Ugur Demirors - 2231819
% Part I, source directly moves towards the sensor.
clear
close all 
N = 6; % Monte carlo trial number
c = 340; % Signal speed.
u = 10/36; % km/h to m/s convertion constant
v = [40 50 60 70 80 90] *u; % Source speed values. (Average car speed values.)
delta  = 400; % delta value ( maximum time value, so max range is 340*2000 = 680km)
f_o = [250 500 750 900 1100 1400]; % Frequency values, baseband signals (low frequency) 
% We inspected very low frequencies because of the performance and ram
% issues.
doppler_multiplier = c ./ (c-v);
Fs =2* f_o .*doppler_multiplier; % Fs equals to Nyquist rate of Doppler shifted f_o.
Ts = 1./Fs;
n = [ 4 5 6 7 8 9]; % n values for t_0
k = [0 1 2 3 4 5]; % k values for t_0
A =5;
m = [ 10 20 30 40 50 60]; % m values. In FM modulation when beta value(proportional with m) increases, signal becomes more 
% noise immuned, but its bandwidth increases (which is bad, it is a trade-off).
first_index = zeros(2*N,1);
last_index = zeros(2*N,1);
real_range = zeros(N,1);
s_r = cell(6,1);
s_r_sin = cell(6,1);
disp('First six signals are linear chirp and last six signals are sin')
for ii=1:N
t = (0:1/Fs(ii):delta-1/Fs(ii)).'; % Time vector
t_0 = 30; % Shifting amount
t_last =360; % range is between t_last and t_0
first_index(ii) = round(t_0*Fs(ii));
last_index(ii) = round(t_last*Fs(ii));
real_range(ii) = (c+v(ii) )* t(first_index(ii)) ; % Determine the real range
noise = 1.2*randn(length(t),1);
s_r{ii} = [zeros(first_index(ii)-1,1);A * cos( doppler_multiplier(ii)* 2 * pi  *( f_o(ii) * t(first_index(ii):last_index(ii)) + (doppler_multiplier(ii)*m(ii) ./ (2 * delta) * (t(first_index(ii):last_index(ii)).^2) ) ) ) ;zeros(length(t)-last_index(ii),1)] + noise;%Linear chirp signal
s_r_sin{ii} = [zeros(first_index(ii)-1,1);A * cos( doppler_multiplier(ii)* 2 * pi  *( f_o(ii) * t(first_index(ii):last_index(ii))  ) ) ;zeros(length(t)-last_index(ii),1)] + noise; % Sin signal
end
all_datas = [s_r; s_r_sin]; % This is for estimation only in single for loop.
Calculation_Fs = [Fs Fs]; % We used same Fs values for sin and linear chirp signals, so Fs values are same.
Calculation_f_o = [f_o f_o];
Calculation_v = [v v];
Calculation_real_range = [ real_range; real_range];
% Part II
PSD_matrix= cell(2*N,1); % To calculate maximum values stft at each frequency
Freq_axis_vector = cell(2*N,1);
Time_axis_vector = cell(2*N,1);
freq_index = zeros(2*N,1);
time_index1 = zeros(2*N,1); % Two indexes for each signal.
estimated_range = zeros(2*N,1);
f_output = zeros(2*N,1);
extracted_speed = zeros(2*N,1);
for ii=1:2*N % We have 2*N signals (sin and linear chirp)
    figure(ii)
    [PSD_matrix{ii},Freq_axis_vector{ii},Time_axis_vector{ii}] = myspectrogram1(all_datas{ii},Calculation_Fs(ii),rectwin(2000),2000,500);
    temp_psd_vector = max(20*log10(abs(PSD_matrix{ii}))); % Find maksimum frequencies at each time instant
    %and determine starting and ending indexes for the signal to  find the range.
    temp_freq_axis_vector = Freq_axis_vector{ii};
    temp_time_vector = Time_axis_vector{ii};
    temp = 0;
    tt = 1;
    while true
        difference = temp_psd_vector(tt+2)- temp_psd_vector(tt);
        if difference >17.2 % We determine starting and ending time indexes according to psd value difference in the vector.
            time_index1(ii) = tt+2;
            break;
        end
        tt = tt+1;
    end
    [~,freq_index(ii)] = max(mean(abs(PSD_matrix{ii}),2));
    f_output(ii) = temp_freq_axis_vector(freq_index(ii));
    extracted_speed(ii) = c*(f_output(ii)-Calculation_f_o(ii))./Calculation_f_o(ii);
    % Error between estimated speed and real speed.
    disp(['Speed estimation error for ', num2str(ii), 'th signal = ', num2str(extracted_speed(ii)-Calculation_v(ii) ),' m/s' ])
    estimated_range(ii) = (c+extracted_speed(ii)) *temp_time_vector(time_index1(ii)); % Range of signals.
    disp(['Range estimation error for ', num2str(ii), 'th signal = ', num2str( Calculation_real_range(ii)-estimated_range(ii) ), ' m'] )
end

% This part is with different window length and overlap length calculations.
figure(13)
myspectrogram1(all_datas{1},Calculation_Fs(1),rectwin(1500),1500,500); % Less length

figure(14)
myspectrogram1(all_datas{1},Calculation_Fs(1),rectwin(1500),1500,500); % More length

figure(15)
myspectrogram1(all_datas{1},Calculation_Fs(1),rectwin(2000),2000,300); % Less overlap length

figure(16)
myspectrogram1(all_datas{1},Calculation_Fs(1),rectwin(2000),2000,800); % More overlap length

