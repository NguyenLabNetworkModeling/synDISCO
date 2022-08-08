function dxdt=EGFR_PYK2_Model(t,x,Qstim,param,inhibitor)

% Q stimulation
EGF=Qstim(1);
HGF=Qstim(2);
caSTAT3=Qstim(3);


% parameters
kc11=param(1);
km11=param(2);
kc12=param(3);
Km12=param(4);
vm21=param(5);
Km21=param(6);
Ki21=param(7);
vm22=param(8);
kc21=param(9);
Km22=param(10);
vs31=param(11);
vm31=param(12);
Km31=param(13);
kd31=param(14);
vm41=param(15);
Km41=param(16);
kc41=param(17);
kc42=param(18);
Km42=param(19);
vm42=param(20);
Km43=param(21);
kd41=param(22);
kc51=param(23);
Km51=param(24);
kc52=param(25);
km52=param(26);
vs61=param(27);
vm61=param(28);
Km61=param(29);
kd61=param(30);
vm71=param(31);
Km71=param(32);
kc71=param(33);
Km72=param(34);
vm72=param(35);
Km73=param(36);
kd71=param(37);
kd72=param(38);
Km74=param(39);
kc81=param(40);
Km81=param(41);
kc82=param(42);
Km82=param(43);
kc91=param(44);
Km91=param(45);
vm91=param(46);
Km92=param(47);
vm11=param(48);
vm51=param(49);
vm81=param(50);
kc101=param(51);
kc102=param(52);
Km101=param(53);
vm101=param(54);
Km102=param(55);
kc43=param(56);
ki41=param(57);
ki101=param(58);
ka121=param(59);
kd121=param(60);
kalEMD1=param(61);
caEGF=param(62);
caHGF=param(63);
EGFRtot=param(64);
STAT3tot=param(65);
Cbltot=param(66);
PTPtot=param(67);
ERKtot=param(68);
ki11=param(69);
ki51=param(70);


Gefit=inhibitor(1); % Gefitinib
PYK2i=inhibitor(2); % PYK2 inhibtor (catalytic inhibitor)
Stattictot=inhibitor(3);% STAT3 inhibitor
EMD=inhibitor(4); % c-MET inhibitor (catalytic inhibitor)

% state variable
pEGFR=x(1); 
EGFRub=x(2); 
PYK2m=x(3); 
PYK2=x(4); 
pPYK2=x(5); 
pSTAT3=x(6);
cMETm=x(7); 
cMET=x(8); 
pcMET=x(9); 
pCbl=x(10);
aPTP=x(11); 
ppERK=x(12);
%notUsed=x(13);
STAT3uStattic=x(14); 
PYK2tot=PYK2+pPYK2;

% algebraic constraints
EGFR = EGFRtot-pEGFR-EGFRub;
STAT3 = STAT3tot-pSTAT3-STAT3uStattic;
Cbl = Cbltot-pCbl;
PTP = PTPtot-aPTP;
ERK = ERKtot-ppERK;
Stattic = Stattictot-STAT3uStattic; 
tactSTAT3=pSTAT3+kalEMD1*caSTAT3;


% Rate equation
d_pEGFR=kc11*(EGF/(1+Gefit/ki11)+caEGF)*EGFR/(km11+EGFR) - (vm11+kc12*aPTP)*pEGFR/(Km12+pEGFR);
d_EGFRub=(vm21+kc21*pCbl)*EGFR/(Km21+EGFR)*Ki21/(Ki21+PYK2tot/(1+PYK2i/ki51)) - vm22*EGFRub/(Km22+EGFRub);
d_PYK2m=vs31 + vm31*tactSTAT3/(Km31+tactSTAT3) - kd31*PYK2m;
d_PYK2=vm41*PYK2m/(Km41+PYK2m) - (kc41*pEGFR+kc42*pcMET/(1+EMD/ki41))*PYK2/(Km42+PYK2) ...
    + (vm42+kc43*aPTP)*pPYK2/(Km43+pPYK2) - kd41*PYK2;
d_pPYK2=(kc41*pEGFR+kc42*pcMET/(1+EMD/ki41))*PYK2/(Km42+PYK2) - (vm42+kc43*aPTP)*pPYK2/(Km43+pPYK2);
d_pSTAT3=kc51*(pPYK2/(1+PYK2i/ki51))*STAT3/(Km51+STAT3) - (vm51+kc52*aPTP)*pSTAT3/(km52+pSTAT3);
d_cMETm=vs61 + vm61*tactSTAT3/(Km61+tactSTAT3) - kd61*cMETm;
d_cMET=vm71*cMETm/(Km71+cMETm) - (kc71*HGF+caHGF)*cMET/(Km72+cMET) + vm72*pcMET/(Km73+pcMET)...
    -(kd71+kd72*pCbl)*cMET/(Km74+cMET);
d_pcMET=(kc71*HGF+caHGF)*cMET/(Km72+cMET) - vm72*pcMET/(Km73+pcMET);
d_pCbl=kc81*pEGFR*Cbl/(Km81+Cbl) - (vm81+kc82*aPTP)*pCbl/(Km82+pCbl);
d_aPTP=kc91*pEGFR*PTP/(Km91+PTP) - vm91*aPTP/(Km92+aPTP);
d_ppERK=(kc101*pcMET/(1+EMD/ki101)+kc102*pEGFR)*ERK/(Km101+ERK)...
    -  vm101*ppERK/(Km102+ppERK);
d_STAT3uStattic =  ka121*STAT3*Stattic - kd121*STAT3uStattic;
notUsed=0;


dxdt=[d_pEGFR;d_EGFRub;d_PYK2m;d_PYK2;d_pPYK2;d_pSTAT3;d_cMETm;d_cMET;...
    d_pcMET;d_pCbl;d_aPTP;d_ppERK;...
    notUsed;d_STAT3uStattic];

end
