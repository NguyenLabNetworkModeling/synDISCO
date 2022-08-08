%% Run 


%% Estimate IC50 value
% load drug_response_matrix_5.mat % (n = 5)
load drug_response_matrix_30.mat % (n = 30)


combo_index = drug_response_matrix.combo_index;
num_of_doses =  drug_response_matrix.num_of_doses;
tspan = drug_response_matrix.tspan;
combo_label = drug_response_matrix.combo_label;
output_labels = drug_response_matrix.output_labels;

for rdout = 1:2 % pERK(1) and pSTAT3(2)
    for ii = 1:6 % all six combinations
        
        IC_D1 = [];
        IC_D2 = [];
        IC_D1(:,1)=linspace(0,drug_response_matrix.nax_conc{rdout}(ii,1),num_of_doses)*I0(combo_index(ii,1));% (EGFRi,pERK)
        IC_D2(:,1)=linspace(0,drug_response_matrix.nax_conc{rdout}(ii,2),num_of_doses)*I0(combo_index(ii,2));% (EGFRi,pERK)
        
        readTime = 24; % hr after treatment
        doseResp(:,:) = drug_response_matrix.drug_response_6d(tspan==readTime*60,rdout,rdout,ii,:,:);
        
        R1(:,1) = doseResp(:,1);
        R2(:,1) = doseResp(1,:);
        
        
        % Calculation of IC50 and plot
        figure('Position',[681   824   508   155])
        
        % Drug 1
        [fitobject1,gof1] = Calculate_IC50_curve(IC_D1,R1);
        IC50(ii,1,rdout) = fitobject1.IC50;
        GOF(ii,1,rdout) = gof1.rsquare;
        
        subplot(1,2,1),plot(fitobject1,IC_D1,R1)
        xlabel(strcat(combo_label{ii,1},'(Conc.)'))
        ylabel(strcat('Activity',{' '},'(',output_labels{rdout},')'))
        grid on
        legend off
        pbaspect([4 3 1])
        
        % Drug 2
        [fitobject2,gof2] = Calculate_IC50_curve(IC_D2,R2);
        IC50(ii,2,rdout) = fitobject2.IC50;
        GOF(ii,1,rdout) = gof2.rsquare;
        
        subplot(1,2,2),plot(fitobject2,IC_D2,R2)
        xlabel(strcat(combo_label{ii,2},'(Conc.)'))
        ylabel(strcat('Activity',{' '},'(',output_labels{rdout},')'))
        grid on
        legend off
        pbaspect([4 3 1])
    end
end

%%

%% SAVE RELEVANT DATA

ic50_matrix.IC50 = IC50;
ic50_matrix.GOF = GOF;
ic50_matrix.combIndex = combo_index;
ic50_matrix.drug_label = combo_label;

fname = 'ic50_matrix.mat';
save(strcat(workdir,'\Outcome\',fname),'ic50_matrix')
