function [Dm,m,gof]=medianEffectEquaiton_Basic(D1,logfa)


%ft = '(x/Dm)^m';
ft = '1/(1+(Dm/x)^m)';


f1=fittype(ft,'coefficients',{'Dm','m'});

% Note: dose = 0 removed
fopt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[100,100],...
    'StartPoint',[3, 1]);
[fo,gof]=fit(D1(2:end),logfa(2:end),f1,fopt);

%plot(fo,log10(D1(2:end)),logfa(2:end))
m=fo.m;
Dm=fo.Dm;

if gof.rsquare < 0.9
    disp('Not well fitted. Try it again with differen option')
    disp(strcat('GOF is ',num2str(gof.rsquare)))
    figure
    plot(fo,D1(2:end),logfa(2:end))
    title(strcat('MEDIAN-EFFECT CURVE / POOR FITTING(R2 = ',num2str(gof.rsquare),')'))
end

