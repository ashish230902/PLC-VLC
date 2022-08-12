%% PLC Model
function stp2=PLC(D,SNR_dB)
%% PLC Model
noise_x =(0:0.001:10);
Pt=1000;
%D = 500;
f=0.02;

%SNR_dB = 5;
% if 75<D<125
%      a0= 9.40*10^(-3);
%      a1=  4.20*10^-7;
%     %b0=(-0.00009098*D)-0.000001126;
%     k=0.7;
% elseif 125<D<175
%     a0= 1.09*10^(-3);
%      a1=  3.36*10^-7;
%     %b0=(-0.00009098*D)-0.000001126;
%     k=0.7;
% elseif 175<D<250
%     a0= 9.33*10^(-3);
%      a1=  3.24*10^-7;
%     %b0=(-0.00009098*D)-0.000001126;
%     k=0.7;
% elseif 250<D<340
%      a0= 8.40*10^(-3);
%      a1=  3*10^-9;
%     %b0=(-0.00009098*D)-0.000001126;
%     k=1;
% elseif 340<D<405
%      a0= 6.20*10^(-3);
%      a1=  4.00*10^-9;
%     %b0=(-0.00009098*D)-0.000001126;
%     k=1;
% end
% 
 a0 = 9.33*10^-3; 
 a1 = 5.1*10^-3; 
 k = 0.7;


h_upperlimit = 100; %Upperlimit of h
h_lowerlimit = 0; %Lowerlimit of h 
x= 0.01; %Interval of PDF
X= (h_lowerlimit:x:h_upperlimit);
X2=X.^2;
A_fn= exp(-2*D*(a0+(a1*(f^k))));
%scale= exp(-1*(a0+(a1*(f.^k)))*D);
scale=sqrt(Pt*A_fn/(pi-2));
%fn = (1/A_fn)*exp(-X2/A_fn);
fn = raylpdf(X,scale);
plot(X2,A_fn)
%% SNR
h2= X.^2;
 Snr_th= 22.6;
 SNR_th=10.^(Snr_th./10);
 %SNR_th= SNR_th.*h2;

 %N0= scale*(sqrt(pi/2));
 N0= (Pt.*h2./SNR_th);
 
%    SNR_dB=30; %SNR value in dB
   SNR=10.^(SNR_dB./10); 
  
  3%SNR= SNR.*h2;
   Noise_var_pdf=sqrt(Pt./(SNR));
   Noise_scale=sqrt(2*Noise_var_pdf/pi);
  
  noise_pdf = raylpdf(noise_x,Noise_scale);
  figure(2)
  plot(noise_x, noise_pdf)
  title('noise pdf')
  area2=trapz(noise_x, noise_pdf)
  
  snr_h= SNR.*h2;
  
  mat= noise_pdf.*(fn.');
  for i=1:length(fn)
  [c limit]=  min(abs(N0(i)-noise_x));
  
       if limit==1
          stp1=0
       else
           
           %stp0= noise_pdf(1:limit).*(GGpdf_h.');
           
            stp1(i)= trapz(noise_x(1:limit),mat(i,1:limit));
            
           %stp= trapz(noise_x(1:limit),noise_pdf(1:limit))*trapz(h,GGpdf_h)
           
       end
  end
  stp2=trapz((X),stp1)

       
 BER=trapz(X,0.5.*erfc(sqrt(snr_h./2)).*fn)

end
