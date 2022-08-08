function CDI = Calculation_CDI_Score(D1,D2,dose_resp)

   for ii = 1:length(D1)
        for jj = 1:length(D2)
            
            % for CDI
            DR1 = dose_resp(ii,1);
            DR2 = dose_resp(1,jj);
            DR12 = dose_resp(ii,jj);
            CDI(ii,jj) = DR12/(DR1*DR2);
          
        end
    end