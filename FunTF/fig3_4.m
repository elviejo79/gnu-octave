time = 0:.001:5;
n = length(time);

frequency = [3 1 6 12];
amplitude = [10 2 5 8];
phase_angle = [0 pi/4 -pi pi/2];

swave = zeros(size(time));
for i=1:length(amplitude)
        swave = swave + amplitude(i)*sin(2*pi*frequency(i)*time+phase_angle(i));
end
        
                                    

clf
plot(time,swave)
xlabel('Time (s)')
ylabel('Amplitude')
  
