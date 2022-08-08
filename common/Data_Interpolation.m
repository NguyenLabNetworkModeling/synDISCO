function [D1_resample,D2_resample,doseRespDat_resample]=Data_Interpolation(D1,D2,doseRespDat,samp_size)

[xData, yData, zData] = prepareSurfaceData( D1, D2, doseRespDat);

% Fit model to data.
[fitresult, ~] = fit( [xData, yData], zData, 'cubicinterp', 'Normalize', 'on' );

D1_resample(:,1)=linspace(min(D1),max(D1),samp_size);
D2_resample(:,1)=linspace(min(D2),max(D2),samp_size);

n1d=length(D1_resample);
n2d=length(D2_resample);

parfor masterIDX=1:(n1d*n2d)
    
    [idx1,idx2]=ind2sub([n1d,n2d],masterIDX);
    doseRespDat_resample_par(masterIDX)=fitresult(D1_resample(idx1),D2_resample(idx2));
    
end

% RESHAPTING
doseRespDat_resample = reshape(doseRespDat_resample_par,[n1d,n2d]);
