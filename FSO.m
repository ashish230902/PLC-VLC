
%% FSO model
%function stp2=FSO(le)
noise_x =(0:0.01:100);
Rdia= 5*10^-3; %aperture of Rx
wav= 1550*10^-9;
h_upperlimit = 10; %Upperlimit of h
h_lowerlimit = 0; %Lowerlimit of h 
x= 0.001; %Interval of PDF
h= (h_lowerlimit:x:h_upperlimit);
Pt_av=100;
k_beam = (2*pi)/wav;
Cn_sq=8.5e-15;
le=3000;
%Cn_sq =6.352*10^-7*le^-2.966;

%  if le<19
%      Cn_sq=0;
%  elseif 19<le<230
%      Cn_sq=4.008*10^(-13)*le^(-1.054);
%  elseif 230<le<850
%      Cn_sq=1.3*10^-15
%  else
%      Cn_sq= 6.352*10^(-16)*le^(-2.966);
%  end

%% Gamma-Gamma

sigma_sq = 1.23*Cn_sq*(k_beam^(7/6))*le^(11/6);

dsig    = sqrt((0.25*k_beam*(Rdia)^2)/le);

%beta1 
beta1_Num = 0.51*sigma_sq*(1+(0.69*(sigma_sq)^1.2))^(-5/6);
beta1_Den = 1+(0.9*dsig^2)+(0.62*dsig^2*(sigma_sq)^1.2);

beta1     = ((exp(beta1_Num/beta1_Den))-1)^-1;

%alpha1
alpha1_Num= 0.49*sigma_sq;
alpha1_Den= (1+(0.65*dsig^2)+(1.11*(sigma_sq)^1.2))^(7/6);

alpha1    = ((exp(alpha1_Num/alpha1_Den))-1)^-1;
GGpdf_h = zeros(1,length(h));

for i = 1:length(h)
    GGpdf_h(i) = PGG2(alpha1,beta1,h(i));
end
% 
% figure(1)
% plot(h,GGpdf_h)
area_pdf=trapz(h, GGpdf_h)

h2=h.^2;

%% SNR

 Snr_th= 22.6;
 SNR_th=10.^(Snr_th./10);
 %SNR_th= SNR_th.*h2;
 n0= sqrt(Pt_av./(SNR_th));
 scale= sqrt(n0/2);
 
 %N0= scale*(sqrt(pi/2));
 N0= (Pt_av.*h2./SNR_th);
 
   SNR_dB=30; %SNR value in dB
   SNR=10.^(SNR_dB./10); 
  
  %SNR= SNR.*h2;
  Noise_var_pdf=sqrt(Pt_av./(SNR));
  Noise_scale=sqrt(Noise_var_pdf/2);
  
  noise_pdf = raylpdf(noise_x,Noise_scale);
%   figure(2)
%   plot(noise_x, noise_pdf)
%   title('noise pdf')
  area2=trapz(noise_x, noise_pdf)
  
  snr_h= SNR.*h2;
  
  mat= noise_pdf.*(GGpdf_h.');
  for i=1:length(GGpdf_h)
  [c limit]=  min(abs(N0(i)-noise_x));
  
       if limit==1
          stp1(i)=0;
       else
         %  stp0= noise_pdf(1:limit).*(GGpdf_h.');
            stp1(i)= trapz(noise_x(1:limit),mat(i,1:limit));
            
           %stp= trapz(noise_x(1:limit),noise_pdf(1:limit))*trapz(h,GGpdf_h)
       end
  end
  stp2=trapz((h),stp1)
       
 BER=trapz(h,0.5.*erfc(sqrt(snr_h./2)).*GGpdf_h)
%end
%%