function [fitobject1,gof1] = Calculate_IC50_curve(D1,R1)

       H0 = 1:2:20;
        % Drug 1
        for kk = 1:length(H0)
            [fitobject1,gof1] = FourParameterLogisticCurve(D1,R1,H0(kk));            
            if gof1.rsquare > 0.9
                break
            end
        end
        
        
        
%     [fitobject1,gof1] = FourParameterLogisticCurve(D1_resample,S1,5);
%     IC50 = fitobject1.IC50;
%     
%     figure('Position',[1189         726         285         258])
%     plot(fitobject1,D1_resample,S1)
%     xline(fitobject1.IC50)
%     yline(1-(1-fitobject1.ymin)/2)
%     xlabel(D1_Name{1})
%     ylabel('Activity')