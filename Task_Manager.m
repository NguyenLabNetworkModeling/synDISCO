%% Task Manager

global param_names state_names inhibitors qstims p0 X0 I0 q0 trainData


% loading data
tbl_param = readtable('BestFittedParam.xlsx',...
    'sheet','parameters',...
    'ReadVariableNames',true,...
    'ReadRowNames',true);
param_names =  tbl_param.name;
p0 =  tbl_param.value;

tbl_x0 = readtable('BestFittedParam.xlsx',...
    'sheet','initial',...
    'ReadVariableNames',true,...
    'ReadRowNames',true);
state_names = tbl_x0.name;
X0 = tbl_x0.value;

tbl_inhibitor = readtable('BestFittedParam.xlsx',...
    'sheet','inhibitor',...
    'ReadVariableNames',true,...
    'ReadRowNames',true);
inhibitors = tbl_inhibitor.name;
I0 = tbl_inhibitor.value;

tbl_qstim = readtable('BestFittedParam.xlsx',...
    'sheet','Qstim',...
    'ReadVariableNames',true,...
    'ReadRowNames',true);
qstims = tbl_qstim.name;
q0 = tbl_qstim.value;

trainData = readtable('modelCalibrationData.xlsx',...
    'ReadVariableNames',true,...
    'ReadRowNames',false);





% IGFR EGFR ERBB PI3K PDK1 IRS
% mTOR
switch jobcode
    case 'model calibration'

        workdir     = strcat(rootwd,'\model calibration');
        mkdir(workdir,'Outcome')
        addpath(genpath(workdir));


        Model_Calibration_Result


    case 'drug response simulation'
        workdir     = strcat(rootwd,'\drug response simulation');
        mkdir(workdir,'Outcome')
        addpath(genpath(workdir));


        % Select a simulation tas from the list
        listJob = {
            'combo inhibition';        % (1)
            'estimate IC50';  % (2)
            };
        [ jobID ] = readInput( listJob );
        sub_jobcode = listJob{jobID}; % Selected

        switch sub_jobcode
            case 'combo inhibition'
                Combined_Inhibition
            case 'estimate IC50'
                Estimate_IC50
        end




    case 'calculate synergy score'
        workdir     = strcat(rootwd,'\calculate synergy score');
        mkdir(workdir,'Outcome')
        addpath(genpath(workdir));


        [file,path,indx] = uigetfile('*.mat','Select [drug_response_matrix]');
        load(file) ;


        listJob = drug_response_matrix.output_labels;
        [ jobID ] = readInput( listJob );
        rd_selectoin = listJob{jobID}; % Selected


        listJob = {
            'scale in a drug concentration';        % (1)
            'scale in a ICx';  % (2)
            };
        icx_selection  = readInput( listJob );


        CalculationDrugSynergy

        
end
