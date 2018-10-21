% Question 2
fs = 44100;       % Sampling frequency (Hz)
t = 0 : 1/fs : 1; % Time axis (seconds)
f = 880;          % Frequency (Hz)
note = {};
vui = [];
for i=1:12
   note{i}=0.02*sin(2*pi*f*(2^(i-1)).^(1/12)*t);
end
do = note{1};
% Do = 'A' 'A#', Re = 'B', Mi = 'C' 'C#', Fa = 'D' 'D#', Son = 'E', La = 'F' 'F#', Si = 'G' 'G#'
A = note{1}; A1 = note{2}; B = note{3}; C = note{4}; C1 = note{5}; D = note{6}; D1 = note{7}; E = note{8}; F = note{9}; F1 = note{10}; G = note{11}; G1 = note{12};
vui = [A A B A D C  A A B A E D ];
[g,fs] = audioread('orig_input.wav');% Reads data from the file named filename, and returns sampled data, y, and a sample rate for that data, Fs.
melody = g + (vui(1:length(g)))'; 
soundsc(melody,44100); % Playback sinewave

% Question 3
Y = fft(g);
plot(abs(Y))

N = length(melody);    % Number of FFT point
transform = fft(g,N)/N;
magTransform = abs(transform);
figure(3);
faxis = linspace(-N/2,N/2,N);
plot(faxis,fftshift(magTransform));
xlabel('Frequency (Hz)')
title('The signal frequency spectrum');
%View frequency content up to half the sampling rate:
  %axis([0 length(faxis)/2, 0 max(magTransform)])
% change the tick labels of the graph from scientific notation to floating point: 
  %xt = get(gca,'XTick');  
  %set(gca,'XTickLabel', sprintf('%.0f|',xt))

%Question 4
figure(4);
win = 128;  % Window length in samples
% Number of samples between overlapping windows:
hop = win/2;          
nfft = win; % width of each frequency bin 
spectrogram(melody,win,hop,nfft,fs,'yaxis')
% Change the tick labels of the graph from scientific notation to floating point: 
yt = get(gca,'YTick');  
set(gca,'YTickLabel', sprintf('%.0f|',yt))
title('Spectrogram of the melody');

