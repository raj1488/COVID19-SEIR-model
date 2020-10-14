clearvars;
close all;
clear all
clc;
mkdir('E:\Desktop content\Publishable manuscript\COVID 19 data\SEIR_model\India data\fitted curves\new_Dir1');
dir1='E:\Desktop content\Publishable manuscript\COVID 19 data\SEIR_model\India data\fitted curves';
dir2='E:\Desktop content\Publishable manuscript\COVID 19 data\SEIR_model';
dir3='E:\Desktop content\Publishable manuscript\COVID 19 data\SEIR_model\India data\fitted curves\new_Dir1';
stateIDVal=[1 23 33 13 11 31 22 36 37 3 6 18 28 30 34 14 16 19 5];
% L=5;
% name3=sprintf('Death_Rate_L%d.xls',L);
% name2=sprintf('Recovery_Rate_L%d.xls',L);
% name1=sprintf('Parameters_Value_L%d.xls',L);
dt = 01; % time step
stDate=datetime(2020,04,15,0,0,0);
endDate=datetime(2020,10,12,0,0,0);
timestart = datetime(2020,03,14,0,0,0):dt:stDate;
start=length(timestart);
StatePopulation1=[417036;53903393;1570458;35607039;124799926;1158473;29436231;615724/2;615724/2;18710922;
                 1586250;63872399;28204692;7451955;13606320;38593948;67562686;35699443;289023;73183;
                 85358965;123144223;3091545;3366710;1239244;2249695;46356334;1413542;30141373;81032689;
                 690251;77841267;39362732;4169794;237882725;11250858;99609303];
StatePopulation(1,1)=sum(StatePopulation1);
StatePopulation(2:length(StatePopulation1)+1,1)=StatePopulation1;

A=csvread('TotalNewCases.csv');
B=csvread('TotalRecoveredCases.csv');
C=csvread('TotalDeathCases.csv');

[~,text1,~]=xlsread('stateName.xlsx');
s=1;
stateName{s}='India';
s=s+1;
stateName{s}=text1{2,2};
for ii=2:length(text1)-1
    if(~strcmp(text1{ii,2},text1{ii+1,2}))
        s=s+1;
        stateName{s}=text1{ii+1,2};
    end
end

for rk=1%:length(stateIDVal)
    rk
    fittedEndDate=datetime(2021,12,07,0,0,0);
time1 = stDate:dt:endDate;
N = numel(time1);
t = [0:N-1].*dt;
stateID=stateIDVal(rk);
Npop=StatePopulation(stateID);
Q=(A(start:N+start-1,stateID)-B(start:N+start-1,stateID)-C(start:N+start-1,stateID))';
R=B(start:N+start-1,stateID)';
D=C(start:N+start-1,stateID)';
Q0 = Q(1); % Initial number of infectious that have bee quanrantined
I0 = 2*Q0; % Initial number of infectious cases non-quarantined
E0 = I0;%A(start); % Initial number of exposed
R0 = B(start,stateID); % Initial number of recovereds
D0 = C(start,stateID); % Initial number of deads
lambdaGuess = [0.01 0.5 0.1];
kappaGuess = [0.01 0.1,10];
alphaGuess = 0.05;
betaGuess = 0.7;
deltaGuess = 0.2;
gammaGuess = 0.3;
guess = [alphaGuess,betaGuess,deltaGuess,gammaGuess,lambdaGuess,kappaGuess]; % initial guess
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,lambdaFun,kappaFun0] = fit_SEIQRDP(Q,R,D,Npop,E0,I0,time1,guess,'Display','off');
R_0(rk)=(beta1/delta1)*((1-alpha1).^(N));
alpha1Val(rk)=alpha1;
betaVal(rk)=beta1;
incubation_period(rk)=1/gamma1;
quarantine_rate(rk)=delta1;
quarantineTime(rk)=1/delta1;
towLambdaVal(rk)=Lambda1(3);
towKappaVal(rk)=Kappa1(3);
% N = numel(time1);
% t = [0:N-1].*dt;
recoveryRate(rk,:)=lambdaFun(Lambda1, t);
deathRate(rk,:)=kappaFun0(Kappa1,t);
dt=01;
time2 = stDate:dt:fittedEndDate;
N = numel(time2);
t = [0:N-1].*dt;
total=Q+R+D;
[S1,E1,I1,Q1,R1,D1,P1] =SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t,lambdaFun,kappaFun0);

figure;

plot(time1,Q,'c',time1,R,'g',time1,D,'r',time1,total,'b','linewidth',2);
hold on
total1=Q1+R1+D1;
plot(time2,Q1,'c-.',time2,R1,'g.',time2,D1,'r--',time2,total1,':b','linewidth',2);
hold on
name1=sprintf(' %s \n Tot Case= %d \n Death= %d \n R0= %.3f',stateName{stateID},max(total1),max(D1), R_0(rk));
 title(name1)
 leg = {'Active','Recovered','Dead','Total Cases','Fitted active cases','Fitted recovered','Fitted Dead','Fitted Total Cases'};
legend(leg{:},'location','best')
 axis tight
% figure;
% plot(time1,Q,'r',time1,R,'c',time1,D,'g',time1,total,'b','linewidth',2);
% hold on
% total1=Q1+R1+D1;
% plot(time1,Q1(1:length(time1)),'k-.',time1,R1(1:length(time1)),'k:',time1,D1(1:length(time1)),'k--',time1,total1(1:length(time1)),'.b','linewidth',2);
% hold on
% name1=sprintf(' %s',stateName{stateID});
%  title(name1)
%  leg = {'Quarantined','Recovered','Dead','Total Cases','Fitted quarantined','Fitted recovered','Fitted Dead','Fitted Total Cases'};
% legend(leg{:},'location','best')
%  axis tight
 cd(dir3)
 name=sprintf('Curve_%s.tiff',stateName{stateID});
 saveas(gcf,name);
%  close all
 cd(dir2)
end





% cd(dir1)
% fid=fopen(name1,'w');
% for ii=1:length(R_0)
%     fprintf(fid,'%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n',R_0(ii),100*alpha1Val(ii),100*betaVal(ii),incubation_period(ii),100*quarantine_rate(ii),quarantineTime(ii),towLambdaVal(ii),towKappaVal(ii));
% end
% fclose(fid);
% fid=fopen(name2,'w');
% for ii=1:length(recoveryRate(:,1))
%     for jj=1:length(recoveryRate(1,:))
%     fprintf(fid,'%.2f\t',100*recoveryRate(ii,jj));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);
% fid=fopen(name3,'w');
% for ii=1:length(deathRate(:,1))
%     for jj=1:length(deathRate(1,:))
%     fprintf(fid,'%.3f\t',100*deathRate(ii,jj));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);