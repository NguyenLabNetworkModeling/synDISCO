%% Run

%% Estimate IC50 value
% load drug_response_matrix_5.mat
% load ic50_matrix_5.mat

load drug_response_matrix_30.mat
load ic50_matrix_30.mat


combo_index = drug_response_matrix.combo_index;
num_of_doses =  drug_response_matrix.num_of_doses;
tspan = drug_response_matrix.tspan;
combo_label = drug_response_matrix.combo_label;
output_labels = drug_response_matrix.output_labels;


rdout  = find(ismember(output_labels,rd_selectoin));
measured_time = 24;

%% USER INPUT

for mm = 1:6
    
    combID = mm; % all six combinations
    
    % drug conc range
    D1{mm}(:,1)=linspace(0,drug_response_matrix.nax_conc{rdout}(combID,1),num_of_doses)*I0(combo_index(combID,1));% (EGFRi,pERK)
    D2{mm}(:,1)=linspace(0,drug_response_matrix.nax_conc{rdout}(combID,2),num_of_doses)*I0(combo_index(combID,2));% (EGFRi,pERK)
    
    ref_ic50(:,:) = ic50_matrix.IC50(combID,:,rdout);
    if icx_selection == 1
        dr1{mm} = D1{mm};
        dr2{mm} = D2{mm};
    elseif icx_selection == 2
        dr1{mm} = D1{mm}/(ref_ic50(1));
        dr2{mm} = D2{mm}/(ref_ic50(2));
    end
    
    
    doseRespDat(:,:) = drug_response_matrix.drug_response_6d(tspan==measured_time*60,rdout,rdout,combID,:,:);
    
    % normalization to basal level (for CDI)
    dose_resp{mm} = doseRespDat/max(doseRespDat(1));
    % drug effect (for BI and CI)
    drug_effect{mm}= 1-dose_resp{mm};
    
    
    CDI{mm} = Calculation_CDI_Score(D1{mm},D2{mm},dose_resp{mm});
    BI{mm} = Calculation_BI_Score(D1{mm},D2{mm},drug_effect{mm});
    BI{mm}(BI{mm} <0) = NaN;
    CI{mm} = Calculation_CI_Score(D1{mm},D2{mm},drug_effect{mm});
    
end

%% Synergy Score at IC50x
IC50 = ic50_matrix.IC50(:,:,rdout);

for mm =1:6
    icx = 1;
    
    D1_IC50 = IC50(mm,1);
    D2_IC50 = IC50(mm,2);
    
    [~, loc_1] = min(abs(D1{mm}-D1_IC50*icx));
    [~, loc_2] = min(abs(D2{mm}-D2_IC50*icx));
    
    SC_Score(mm,1) = log2(CDI{mm}(loc_1,loc_2)); % CDI
    SC_Score(mm,2) = log2(BI{mm}(loc_1,loc_2)); % BI
    SC_Score(mm,3) = log2(CI{mm}(loc_1,loc_2)); % CI
    
end



fig_1 = figure('position',[681   721   368   258]);
% Score was normalized to its maxmum
bar(SC_Score./repmat(max(abs(SC_Score)),size(SC_Score,1),1))
legend({'CDI','BI','CI'})
xticklabels(combo_label)
xtickangle(45)
ylabel(strcat('Synergy Score',{' '},'(',output_labels{rdout},')'))
title(strcat('ICx = ',{''},num2str(icx)))


fname_fig_1 = strcat(workdir,'\Outcome','\','SynergyScore_',...
    output_labels{rdout},'_',...
    'IC_',num2str(icx),'.jpeg');
saveas(fig_1,fname_fig_1)


for ii = 1:size(combo_label,1)
    drug_pairs{ii,1} = strcat(combo_label{ii,1},'+',combo_label{ii,2});
end

% Save results
SC_score_SIM = array2table(SC_Score,'VariableNames',{'CDI(Log2)','BI(Log2)','CI(Log2)'},...
    'RowNames',drug_pairs);



fname = strcat(workdir,'\Outcome','\','SynerScore_SIM_IC_',num2str(rdout),'_',...
    num2str(measured_time),'_',num2str(icx),'.xlsx');
writetable(SC_score_SIM,fname,'WriteVariableNames',true,...
    'WriteRowNames',true)



%% PLOTS

for mm = 1:6 % six combinations
    
    plt_position = [18   343   570   476];
    fig_2 = figure('position',plt_position);
    sc = surfc(dr1{mm},dr2{mm},dose_resp{mm}');
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    zlabel(strcat('drugResponse',{' '},'(',output_labels{rdout},')'))
    view([144.4947 56.1727])
    colorbar
    
    fname_fig_2 = strcat(workdir,'\Outcome','\','Drug Response_',...
        output_labels{rdout},'_',...
        string(combo_label(mm,1)),'_',...
        string(combo_label(mm,2)),'_',...
        'IC_',num2str(icx),'.jpeg');    
    saveas(fig_2,fname_fig_2)    
    
    
    
    fig_3 = figure('position',[596         342        1260         471]);
    subplot(2,3,1),
    surfc(dr1{mm},dr2{mm},log2(CDI{mm})')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    zlabel(strcat('Log2CDI',{' '},'(',output_labels{rdout},')'))
    view([147.265179814664 24.3031140777964]);
    colorbar
    
    
    subplot(2,3,2),surfc(dr1{mm},dr2{mm},log2(BI{mm})')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    zlabel(strcat('Log2BI',{' '},'(',output_labels{rdout},')'))
    view([147.265179814664 24.3031140777964]);
    colorbar
    
    subplot(2,3,3),surfc(dr1{mm},dr2{mm},log2(CI{mm})')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    zlabel(strcat('Log2CI',{' '},'(',output_labels{rdout},')'))
    view([147.265179814664 24.3031140777964]);
    colorbar
    
    
    
    subplot(2,3,4),contourf(dr1{mm},dr2{mm},log2(CDI{mm})','ShowText','on')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    title('CDI(Log2)')
    pbaspect([4 3 1])
    colorbar
    
    subplot(2,3,5),contourf(dr1{mm},dr2{mm},log2(BI{mm})','ShowText','on')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    title('BI(Log2)')
    pbaspect([4 3 1])
    colorbar
    
    subplot(2,3,6),contourf(dr1{mm},dr2{mm},log2(CI{mm})','ShowText','on')
    xlabel(string(combo_label(mm,1)))
    ylabel(string(combo_label(mm,2)))
    title('CI(Log2)')
    pbaspect([4 3 1])
    colorbar
    
    fname_fig_3 = strcat(workdir,'\Outcome','\','Drug Synergy_',...
        output_labels{rdout},'_',...
        string(combo_label(mm,1)),'_',...
        string(combo_label(mm,2)),'.jpeg');
    
    saveas(fig_3,fname_fig_3)
    
end





