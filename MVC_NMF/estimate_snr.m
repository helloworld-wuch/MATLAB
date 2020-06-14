function snr_est = estimate_snr(R,r_m,x)  %²âÊÔÔëÉù±È

         [L , ~]=size(R);           % L number of bands (channels)
                                  % N number of pixels (Lines x Columns) 
         [p, N]=size(x);           % p number of endmembers (reduced dimension)

         P_y = sum(R(:).^2)/N;
         P_x = sum(x(:).^2)/N + r_m'*r_m;
         snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
return;

