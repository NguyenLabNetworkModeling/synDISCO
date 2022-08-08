function BI = Calculation_BI_Score(D1,D2,drug_effect)

   for ii = 1:length(D1)
        for jj = 1:length(D2)
            
            % for BI
            DE1 = drug_effect(ii,1);
            DE2 = drug_effect(1,jj);
            DE12 = drug_effect(ii,jj);
            
            BI_prd = DE1 + DE2 - DE1*DE2;
            BI(ii,jj) = BI_prd/DE12;
        end
    end