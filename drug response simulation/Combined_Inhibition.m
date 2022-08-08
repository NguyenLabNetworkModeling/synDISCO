%% Run 

%%
% COMINH_tspan =linspace(0,24,300)*60;

odeopts=[];
% inhibitorLabels=INHIBITOR_val.Properties.RowNames;
num_of_doses= 30; 
% 30 the number of sample size
% 5

%% initialization of ODE
tspan0 = linspace(0,1e6,100);
[~,y0]=ode15s(@EGFR_PYK2_Model,tspan0,X0,odeopts,...
    q0*0,p0,I0*0);


%% Consider all combinations (6 pairs)
% 1: EGFRi(Gef),
% 2: PYK2(PF396),
% 3: Stat2i(STattic),
% 4:cMETi(EMD)

combo_index=[1 2;        % 1: EGFRi(Gef),        2: PYK2(PF396)
    1 3;                % 1: EGFRi(Gef),        3: Stat2i(STattic)
    1 4;                % 1: EGFRi(Gef),        4:cMETi(EMD)
    2 3;                % 2: PYK2(PF396),       3: Stat2i(STattic),
    2 4;                % 2: PYK2(PF396),       4:cMETi(EMD),
    3 4];               % 3: Stat2i(STattic),   4:cMETi(EMD)

combo_label = {'EGFRi', 'PYK2i'
    'EGFRi', 'STAT3i'
    'EGFRi', 'cMETi'
    'PYK2i', 'STAT3i'
    'PYK2i', 'cMETi'
    'STAT3i', 'cMETi'};


% Maximum drug concentration
nax_conc{1} = [0.3 0.25;    % EGFRi PYK2i
    0.3 0.25;                   % EGFRi STAT3i
    0.3 0.05;                   % EGFRi cMETi
    0.3 0.3;                    % PYK2i STAT3i
    0.25 0.06;                  % PYK2i cMETi
    0.25 0.05];                 % STAT3i cMETi

nax_conc{2} = [0.08 0.06;   % EGFRi PYK2i
    0.08 0.15;                  % EGFRi STAT3i
    0.08 0.01;                  % EGFRi cMETi
    0.06 0.15;                  % PYK2i STAT3i
    0.06 0.01;                  % PYK2i cMETi
    0.15 0.01];                 % STAT3i cMETi

output_labels = {'pERK', 'pSTAT3'};

%% Simulation of Drug response

n1d=2; % readout (pERK, pSTAT3) kk --> idx1
n2d=length(combo_index(:,1)); % Num of Combinations mm --> idx2
n3d=num_of_doses; % Drug 1, jj --> idx3
n4d=num_of_doses; % Drug 2, ii --> idx4

tspan =[0 0.5 1:48]*60;

parfor masterIDX=1:(n1d*n2d*n3d*n4d)
    
    disp(masterIDX)
    
    % index to sub
    [idx1,idx2,idx3,idx4]=ind2sub([n1d,n2d,n3d,n4d],masterIDX);
    
    CB_D1=linspace(0,nax_conc{idx1}(idx2,1),num_of_doses)*I0(combo_index(idx2,1));% (EGFRi,pERK)
    CB_D2=linspace(0,nax_conc{idx1}(idx2,2),num_of_doses)*I0(combo_index(idx2,2));% (EGFRi,pERK)
    
    blocker=zeros(1,4);
    % blocker=[EGFRi PYK2i STAT3i cMETi]
    blocker(combo_index(idx2,1))=CB_D1(idx3); % frist drug
    blocker(combo_index(idx2,2))=CB_D2(idx4); % second drug
    
    
    tspan0 = tspan; % simulaiton time = 24 hours
    EGF = q0(1);
    [t,y1]=ode15s(@EGFR_PYK2_Model,tspan0,y0(end,:),odeopts,....
        [EGF 0 0],p0,blocker);
    
    readouts(:,:,masterIDX) = [y1(:,ismember(state_names,'ppERK')) ...
        y1(:,ismember(state_names,'pSTAT3'))];
end

% RESHAPTING THE OUTCOME
[nrow,ncol,~] = size(readouts);
drug_response_6d = reshape(readouts,[nrow,ncol,n1d,n2d,n3d,n4d]);
% {time, variables(pERK,pSTAT3), markers(pERK,pSTAT3), combinations,
% drug1, drug2};

%%

% SAVE RELEVANT DATA
drug_response_matrix.drug_response_6d = drug_response_6d;
drug_response_matrix.nax_conc = nax_conc;
drug_response_matrix.tspan = tspan;
drug_response_matrix.combo_index = combo_index;
drug_response_matrix.combo_label = combo_label;
drug_response_matrix.output_labels = output_labels;
drug_response_matrix.note = ...
    {'time (50)','variable(pERK, pSTAT3)','readout(pERK, pSTAT3)','combinations (6 x 2)','drug1 (30)','drug2 (30)'};
drug_response_matrix.num_of_doses = num_of_doses;

fname = 'drug_response_matrix.mat';
save(strcat(workdir,'\Outcome\',fname),'drug_response_matrix')

%% Plot
% e.g., pERK(1) as readout and EGFRi(Gef) and PYK2i(PF396)combination

readTime = 24;

for ii = 1:2 % {1: pERK, 2: pSTAT3}
    for jj = 1:6 % six combinations
        
        output_sel = ii; % {1: pERK, 2: pSTAT3}
        drugComb = jj;
        
        dat1 = [];
        dat1(:,:) = drug_response_6d((tspan==readTime*60),output_sel,output_sel,drugComb,:,:);
        CB_D1 = linspace(0,nax_conc{output_sel}(drugComb,1),num_of_doses)*I0(combo_index(drugComb,1));% (EGFRi,pERK)
        CB_D2 = linspace(0,nax_conc{output_sel}(drugComb,2),num_of_doses)*I0(combo_index(drugComb,2));% (EGFRi,pERK)
        
        % 3D DOSE-RESPONSE CURVE
        figure('position',[ 218         677        1006         246])
        subplot(1,3,1),surfc(CB_D1,CB_D2,dat1')
        xlabel(string(combo_label(drugComb,1)))
        ylabel(string(combo_label(drugComb,2)))
        zlabel(strcat(output_labels{output_sel},{' '},'(t=',num2str(readTime),'hr)'))
        view([134.4 28.2])
        %[aa,bb]=view
        
        % 2D DOSE-DEPENDENT TIME PROFILES (drug 1)
        dat2 = [];
        dat2(:,:) = drug_response_6d(:,output_sel,output_sel,drugComb,:,1);
        
        subplot(1,3,2),
        plot(tspan/60,dat2')
        xlabel('time (hr)')
        ylabel(output_labels{output_sel})
        title(strcat('Varying ',{' '},combo_label(drugComb,1)))
        box off
        pbaspect([4 3 1])
        xline(24,'--r',{'Measurmet','time'})
        
        % 2D DOSE-DEPENDENT TIME PROFILES (drug 2)
        dat2 = [];
        dat2(:,:) = drug_response_6d(:,output_sel,output_sel,drugComb,1,:);
        
        subplot(1,3,3),
        plot(tspan/60,dat2')
        xlabel('time (hr)')
        ylabel(output_labels{output_sel})
        title(strcat('Varying ',{' '} ,combo_label(drugComb,2)))
        box off
        pbaspect([4 3 1])
        xline(24,'--r',{'Measurmet','time'})
        
    end
end

