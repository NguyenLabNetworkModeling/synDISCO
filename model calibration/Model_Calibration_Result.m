
%% initialization of ODEs
odeopts=[];
CAL_tspan0 = linspace(0,1e6,100);

[~,y0]=ode15s(@EGFR_PYK2_Model,CAL_tspan0,X0,odeopts,...
    q0*0,p0,I0*0);


%% EGF stimulation
tspan = trainData.Time * 60;

[~,y1]=ode15s(@EGFR_PYK2_Model,tspan,y0(end,:),odeopts,...
    q0,p0,I0*0);

out.val = [y1(:,ismember(state_names,'pEGFR')) ...
    y1(:,ismember(state_names,'pSTAT3')) ...
    sum(y1(:,ismember(state_names,{'cMET','pcMET'})),2) ...
    sum(y1(:,ismember(state_names,{'PYK2','pPYK2'})),2) ...
    y1(:,ismember(state_names,'ppERK'))];
out.label = {'pEGFR','pSTAT3','cMET','PYK2','pERK'};

%% Produce the plot
fig_1 = figure('position', [483   659   499   239]);
for ii = 1:length(out.label)
    exp_dat = table2array(trainData(:,ismember(trainData.Properties.VariableNames',out.label{ii})));
    exp_dat = exp_dat/max(exp_dat); % normalied to max
    
    sim_dat = out.val(:,ii);
    sim_dat = sim_dat/max(sim_dat); % normalized to max
    
    subplot(2,3,ii),
    p1 = plot(tspan/60,[exp_dat sim_dat],'-o','LineWidth',1,'MarkerSize',6);
    set(p1(1),'MarkerFaceColor',[0    0.4470    0.7410]);
    set(p1(2),'MarkerFaceColor',[0.8500    0.3250    0.0980]);
    set(gca,'LineWidth',1)
    xlabel('Time (hour)'),
    ylabel(strcat(out.label{ii},{' '},'(normal)')),
    box off
    pbaspect([4 3 1])
    if ii==5
        legend('Exp','Sim')
    end
    
end

fname_fig_1 = strcat(workdir,'\Outcome','\','Model Calib Result.jpeg');
saveas(fig_1,fname_fig_1)


