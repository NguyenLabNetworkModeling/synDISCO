function [Dm,m,gof]=medianEffectEquaiton(D1,logfa)    

ft = fittype('m*x+b');
[fo,gof]=fit(log10(D1(2:end)),logfa(2:end),ft,'StartPoint',[3,-3]);
%plot(fo,log10(D1(2:end)),logfa(2:end))
m=fo.m;
Dm=10^(-fo.b/m);
if gof.rsquare < 0.9
    disp('Not well fitted. Try it again with differen option')
    disp(strcat('GOF is ',num2str(gof.rsquare)))
    figure
    plot(fo,log10(D1(2:end)),logfa(2:end))
    title(strcat('MEDIAN-EFFECT CURVE / POOR FITTING(R2 = ',num2str(gof.rsquare),')'))
    
end

