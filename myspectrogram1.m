function [S1,f,t]= myspectrogram1(signal, fs, window, length_of_window,overlap)
%win=@window;
[S1, f, t]=stft(signal,fs,'Window',window,'OverlapLength',overlap,'FFTLength',length_of_window);
Row=ceil(numel(f)/2);
f=f(Row:end,1);%Take positive frequencies like spectrogram.
S1=S1(Row:end ,:);
surf(t, f, 20*log10(abs(S1)))
shading interp
axis tight
view(0, 90)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Amplitude spectrogram of the signal')
hcol = colorbar;
set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(hcol, 'Magnitude, dB')
end