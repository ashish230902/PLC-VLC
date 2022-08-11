%% This programe provides SER analysis of amplify and forward cooperative protocol using QPSK modulation%%%
%% Results are compared with theoratical results for QPSK modulation with AWGN and Rayleigh fading channels %%%
clear all
clc
N= 10^6; %number of symbols
Pt=1;% Total Tranmited Power
Pb=Pt/2;% Transmitted power of source
Pr=Pt/2; % Transmitted power of Relay
snr_db=[0:40];% SNR in dB
data_BS1 = ((2*(rand(1,N)>0.5)-1) + j*(2*(rand(1,N)>0.5)-1));%% Data generation
data=zeros(1,N); % Received data matrix

%% Channel gains
        channel_BS1_RS1=(1/sqrt(2)).*[randn(1,N)+j*randn(1,N)];% Channel coefficients for Soure-Relay Link
        a11=sum((abs(channel_BS1_RS1)))/N; %Average channel gain for Soure-Relay Link Channel
                
        channel_RS1_MS1=(1/sqrt(2)).*[randn(1,N)+j*randn(1,N)];% Channel coefficients for Relay-Desination Link
        b11=sum((abs(channel_RS1_MS1)))/N; %Average channel gain for Relay-Desination Link Channel
                
        channel_BS1_MS1=(1/sqrt(2))*[randn(1,N)+j*randn(1,N)];% Channel coefficients for Soure-Desination Direct Link
        c11=sum((abs(channel_BS1_MS1)))/N; %Average channel gain for Soure-Desination Link Channel
 %% AWGN     
        nbr=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN for at Relay for Soure-Relay Link
        nrm=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN for at Desination  for Relay-Desination Link
        nbm=(1/sqrt(2)).*(randn(1,N)+j*randn(1,N));%AWGN for at Desination  for source-Desination Link
        
for       i=1:length(snr_db)
                snr_linear(i)=10.^(snr_db(i)/10);   % Linear value of SNR                           
                No(i)=Pt./snr_linear(i);            % Noise power
                 
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr./(a11^2.*Pb+No(i)));
                              
   %% RECEIVED DATA AT RELAY FROM BS         
               data_BS1_RS1=sqrt(Pb/2).*data_BS1.*channel_BS1_RS1+sqrt(No(i))*nbr;
                            
   %% RECEIVED DATA AT MS FROM RS     
        data_RS1_MS1=Beta1.*data_BS1_RS1.*channel_RS1_MS1+sqrt(No(i))*nrm;%+(5.7712e-014 -1.8746e-014i);
        
   %% RECEIVED DATA AT MS FROM BS          
        data_BS1_MS1=sqrt(Pb/2).*data_BS1.*channel_BS1_MS1+sqrt(No(i))*nbm;%+( 5.8496e-008 +4.9957e-007i);
        
   %%   MRC AT RECEIVER 
        data_rd=data_RS1_MS1.*conj(channel_BS1_RS1.*channel_RS1_MS1)+ data_BS1_MS1.*conj(channel_BS1_MS1);
       
   %% DEMODULATION
        g=data_rd;
        c=real(g);
        d=imag(g);
        data(find(c>=0 & d>=0))=1+1*j;
        data(find(c>=0 & d<0))=1-1*j;
        data(find(c<0 & d>=0))=-1+1*j;
        data(find(c<0 & d<0))=-1-1*j;     
        error_af(i)=size(find((data_BS1- data)),2); %% CALCULATING ERRORS
    
end 
%% PLOTTING THE SIMULATION AND THEORATICAL RESULTS
    figure
    %% SIMULATION RESULTS
    simber_af=error_af/N;
    semilogy(snr_db,simber_af,'g.-')
    hold on
    %% THEORATICAL RESULTS (Ref: Digital communication Over Fading Channels  By Alouini) 
    theorySer_fad=(3/4)*[1-sqrt(0.5.*snr_linear./(1+0.5.*snr_linear)).*4/(3*pi).*((pi/2)+atan(sqrt(0.5*snr_linear./(1+0.5*snr_linear))))];
    semilogy(snr_db,theorySer_fad,'r.-'); % Theoretical QPSK SER for fading channel
    hold on
    theorySer_awgn = erfc(sqrt(0.5*(10.^(snr_db/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(snr_db/10))))).^2;
    semilogy(snr_db,theorySer_awgn,'b.-'); % Theoretical QPSK SER for AWGN channel
    %% Theoretical rsult for Amplify and forward (Ref: Cooperative networking ny K.J Ray)
    a1=(var(channel_BS1_RS1));% variance of Soure-Relay Link Channel
    b1=(var(channel_RS1_MS1));% variance of Relay-Destination Link Channel
    c1=(var(channel_BS1_MS1));% variance of Soure-Destination Link Channel
    p_th=(( 0.3608*No.^2)/0.25).*(1/Pb).*((1/Pb)+(1/Pr));
    semilogy(snr_db,p_th,'-*');
    grid on
    axis([0 40 10^-5 0.5]);
    xlabel('SNR(dB)');
    ylabel('SER');
    title('Amplify and Forward Cooperative Communication')
    legend('Simulation-fad','Theory-fad','Theory-awgn','Theory-amp')
    
    
