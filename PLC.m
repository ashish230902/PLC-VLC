%% PLC Model
function stp=PLC(D)
noise_x =(0:0.01:60);

Pt=1000;

% a0= 0 ;
% a1=  7.8*10^-11;
%b0=(-0.00009098*D)-0.000001126;
%k=1;

a0 = (0.0002086*D) +0.0008739; 
a1 = (0.00002644*D) +0.00004644; 
k = (-0.00009098*D)-0.000001126;

f=20;

h_upperlimit = 60; %Upperlimit of h
h_lowerlimit = 0; %Lowerlimit of h 
x= 0.01; %Interval of PDF
X= [h_lowerlimit:x:h_upperlimit];

A_fn= exp(-D*(a0+(a1*(f^k))));
%scale= exp(-1*(a0+(a1*(f.^k)))*D);
scale=sqrt(Pt*A_fn/(pi-2));
fn = raylpdf(X,scale);

plot(X,fn)

%% SNR
h2= X.^2;
 Snr_th= 22.6;
 SNR_th=10.^(Snr_th./10);
 %SNR_th= SNR_th.*h2;

 %N0= scale*(sqrt(pi/2));
 N0= (Pt.*h2./SNR_th);
 
   SNR_dB=30; %SNR value in dB
   SNR=10.^(SNR_dB./10); 
  
  3%SNR= SNR.*h2;
   Noise_var_pdf=sqrt(Pt./(SNR));
   Noise_scale=sqrt(Noise_var_pdf/2);
  
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
