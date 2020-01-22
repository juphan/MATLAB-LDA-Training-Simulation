function y = wav_calc(emg)

dx = emg(1:99) - emg(2:100);
y = mean(abs(dx));

end