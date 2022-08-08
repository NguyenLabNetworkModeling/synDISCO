function [CI,gof1,gof2] = Calculation_CI_Score(D1,D2,drug_effect)
% CI (Pharmacol Rev 58:621–681, 2006)
R1(:,1) = drug_effect(:,1);
R2(:,1) = drug_effect(1,:);

% fraction affected (fa)
logfa_1=log10(R1./(1-R1));
logfa_2=log10(R2./(1-R2));


if isempty(find(imag(logfa_1)))
    [Dm1,m1,gof1]=medianEffectEquaiton(D1,logfa_1);
else
    [Dm1,m1,gof1]=medianEffectEquaiton_Basic(D1,R1);
end

if isempty(find(imag(logfa_2)))
    [Dm2,m2,gof2]=medianEffectEquaiton(D2,logfa_2);
else
    [Dm2,m2,gof2]=medianEffectEquaiton_Basic(D2,R2);
end


% [Dm1,m1]=medianEffectEquaiton(D1,logfa_1);
% [Dm2,m2]=medianEffectEquaiton(D2,logfa_2);

for ii = 1:size(drug_effect,1)
    for jj = 1:size(drug_effect,1)
        % drug effect
        %dxx = 0.9;
        dxx = drug_effect(ii,jj);
        
        % Combined Drug effect
        [~, loc_3] = min(abs(drug_effect(:)-dxx));
        [idx1,idx2]=ind2sub(size(drug_effect),loc_3);
        
        Dx_d1=Dm1*(dxx/(1-dxx))^(1/m1);
        Dx_d2=Dm2*(dxx/(1-dxx))^(1/m2);
        CI(ii,jj)=D1(idx1)/Dx_d1 + D2(idx2)/Dx_d2;
        if dxx < 0.01
            CI(ii,jj)= NaN;
        end
    end
end