function GGpdf = PGG2(alpha,beta,h)

    m1        = round(alpha - beta);
    m2        = ((alpha + beta)/2)-1;
    
    if h == 0
        
        GGpdf = 0;
        disp('----------------');
        disp('h cannot be zero');
        disp('----------------');
    else
        
        GGpdf_Num   = (2*(alpha*beta)^(m2+1))*(h^m2)*besselk(m1,(2*sqrt(alpha*beta*h)));
        GGpdf_Den   = gamma(alpha)*gamma(beta);
        GGpdf       = GGpdf_Num/GGpdf_Den;
        
    end
    
end
