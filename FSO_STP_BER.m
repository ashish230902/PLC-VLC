clear
clc

%*************Simulation Parameters***************************

Nsc= 128; % Number of OFDM subcarrier 
Nsd= 110; % number of used sub carrier
Nsp= 18; %number of used pilot sub-carrier
Rcp=1/9; % Overhead of cyclic prefix 
las_width= 5e6; %width of laser OFDM
Rb= 10e9; %Data Rate
Rin= -150; %Relative Intensity noise db/hz
wav= 1550e-9; %Wavelength 
G= 10; %preamplifying gain db;
n= 0.75; % Effeciency of Photo Diode
Fn= 4.77; %Noise Figure of Optical Amplifier
m = 0.1; %modulation technique
length = 1000; %link length
Tn = 0.9; % Transmitter optical effeciency
Rn= 0.75; %Reciever optical effeienyy 
Tdia = 0.1; %aperture of Tx
Rdia= 0.3; %aperture of Rx
theta= 72e-6; % Half Beam divergence angle
coeff= 2000; % Attenuation coeffecient for haze atmospheric condition
loss_window= 1; %optical Window loss
loss_pointing= 1; % Pointing loss
Ref_in_mod= 1.7e-14;% Refractive index for moderate turbulance
Ref_in_strong= 5e-14;% Refractive index for strong turbulance
G=8/theta;

Pt_av=;


%*************OFDM Simulation***************************
M = 4;                          %   QPSK signal constellation
no_of_data_points = 128;        %   have 64 data points
block_size = 18;                 %   size of each ofdm block
cp_len = 2;  %   length of cyclic prefix
no_of_ifft_points = block_size;           %   8 points for the FFT/IFFT
no_of_fft_points = block_size;
%   ---------------------------------------------
%   B:  %   +++++   TRANSMITTER    +++++
%   ---------------------------------------------
%   1.  Generate 1 x 64 vector of data points phase representations
data_source = randsrc(1, no_of_data_points, 0:M-1);
figure(1)
stem(data_source); grid on; xlabel('data points'); ylabel('transmitted data phase representation')
title('Transmitted Data "O"')
%   2.  Perform QPSK modulation
qpsk_modulated_data = pskmod(data_source, M);
scatterplot(qpsk_modulated_data);title('qpsk modulated transmitted data')
%   3.  Do IFFT on each block
%   Make the serial stream a matrix where each column represents a pre-OFDM
%   block (w/o cyclic prefixing)
%   First: Find out the number of colums that will exist after reshaping
num_cols=length(qpsk_modulated_data)/block_size;
data_matrix = reshape(qpsk_modulated_data, block_size, num_cols);
%   Second: Create empty matix to put the IFFT'd data
cp_start = block_size-cp_len;
cp_end = block_size;
%   Third: Operate columnwise & do CP
for i=1:num_cols,
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    %   Compute and append Cyclic Prefix
    for j=1:cp_len,
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    %   Append the CP to the existing block to create the actual OFDM block
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end
%   4.  Convert to serial stream for transmission
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal');grid on;
P=Pt_av*Tn*G*(1+(m*ofdm_signal);







%   ------------------------------------------
%   E:  %   +++++   RECEIVER    +++++
%   ------------------------------------------
%   1.  Pass the ofdm signal through the channel
recvd_signal = ofdm_signal;
%   4.  Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);
%   5.  Remove CP
recvd_signal_matrix(1:cp_len,:)=[];
%   6.  Perform FFT
for i=1:cols_ifft_data,
    %   FFT
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);
end
%   7.  Convert to serial stream
recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));
%   8.  Demodulate the data
qpsk_demodulated_data = pskdemod(recvd_serial_data,M);
scatterplot(qpsk_modulated_data);title('qpsk modulated received data')
figure(5)
stem(qpsk_demodulated_data,'rx');
grid on;xlabel('data points');ylabel('received data phase representation');title('Received Data "X"')  



