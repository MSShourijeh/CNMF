% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % This is the main file that constructs the matrices for (C)NMF % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function Main_CNMF
close all; 
clear all; clc;

address=pwd;

task='RA';
trial_prefix='normalised_data_stride_';

NoM=12;             % # of muscles
NoD=100;            % # of data points
NoS=4;              % # of synergies
sub_range=1:6;      %
Nsub=1;             % # of subjects
Ntri=numel(sub_range)*Nsub; % Total # of trials
samp_rate=1500;

plot_sub_no=6;      % # of subjects to be plotted together in the reconstruction figure

if NoS==4
    wsize=1;
else
    wsize=NoS/4;
end

load('Sub1_RA_trials_stride.mat')
data_all=zeros(NoD*Ntri,12);
for i=1:6
   temp_data=eval([trial_prefix, num2str(i),'1']);
   t=linspace(1,NoD,length(temp_data));t=t';
   temp_data_interp=interp1(t,temp_data,1:NoD);
   data_all(1+NoD*(i-1):NoD*i,:)=temp_data_interp;
end


%% DATA
%% Zeroing the -ve Values
[aaa,bbb]=find(data_all<0);
for iii=1:length(aaa)
data_all(aaa(iii),bbb(iii))=0;
end


%% Call the CNMF function
for rep=1:1 % number of replicates
[Ctot,S,D,tot_VAF,tot_VAF_all,mean_tot_VAF_all,tot_VAF_Torres] = CNMF(Ntri,data_all,NoS,NoD,NoM);
end
S_final=S;

%% Plot the results
contr = PLOT_CNMF(wsize,Ntri,NoS,S,NoM,plot_sub_no,NoD,Ctot,data_all);
%% Total Contribution of each Synergy
for i=1:NoS
tot_contr(i)=sum(sum(contr(:,:,i)))/sum(sum(sum(contr(:,:,:))));
end
%% ICC for subject comparison: There will be NoS number of ICCs at the end.
type='C-k'; % 2-way mixed, consistency (no interaction)
alpha=0.05;

addpath([pwd '/ICC/']);

for syn=1:NoS
    C_ICC_temp=Ctot(:,syn);
    for sub=1:Ntri
        C_ICC(:,sub)=C_ICC_temp((sub-1)*NoD+1:sub*NoD,:);
    end
    
[r1(syn), LB1(syn), UB1(syn), F1(syn), df11(syn), df12(syn), p1(syn)] = ICC(C_ICC, type, alpha);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ctot,S,CF1,tot_VAF,tot_VAF_all,mean_tot_VAF_all,tot_VAF_Torres]  =  CNMF(Ntri,data_all,NoS,NoD,NoM)
% 5 steps are done in the CNMF function:
% 0) random inputs
% 1) MATLAB nnmf command
% 2) MATLAB lsqnonneg
% 3) MATLAB fmincon
% 4) purturb the results

% [options, options2]  =  optim_Settings; % loads settings for optimization

%% A  =  CS

N_rep  =  3; norm_l  =  'fro';%2;
for i  =  1:N_rep
    
    noc  =  NoS; maxiter  =  500; NoC = NoM;NoR = NoD*Ntri;
    opt  =  statset('MaxIter',maxiter,'Display','final','TolFun',1e-9,'TolX',1e-9);
    if i == 1
        % 0) random inputs
        C = rand(NoR,noc);S = rand(noc,NoC);
    else
        % 4) purturbation
        C = C+rand(size(C));S = S+rand(size(S));
    end
    % 1) nnmf
    [Ctot,temp_S,CF1]  =  nnmf(data_all,noc,'rep',10,...
        'options',opt,...
        'w0',abs(C),'h0',S,'alg','als'); 
    
%     for k = 1:NoM
%         if k =  = 1
%             S = zeros(NoS,NoM);
%         end
%         % 2) lsqnonneg
%         [S_temp,resnorm,residual,exitflag,output,lambda]  =  lsqnonneg(Ctot,data_all(1:NoD*Ntri,k),options2);
%         CF2 = resnorm;
%         S(:,k) = S_temp;
%     end
    
    % 3) fmincon
    lb2 = zeros(NoD*Ntri,NoS);
%     [Ctot,CF1,exitflag,output,lambda,grad,hessian]  =  fmincon(@CF_CNMF_fmincon,Ctot,[],[],[],[],lb2,[],[],options2);
    % 3') fminsearch
%     [Ctot,CF1,exitflag,output] = fminsearchbnd(@CF_CNMF_fmincon,Ctot,lb2,[],options);
    
    S = temp_S;
    C = Ctot;
    
end

%% VAF Definitions
for i = 1:Ntri
    tot_VAF(i) = (1-(var(data_all(NoD*(i-1)+1:NoD*i,:)-Ctot(NoD*(i-1)+1:NoD*i,:)*S)/var(data_all(NoD*(i-1)+1:NoD*i,:))))*100;
end
tot_VAF_all = (1-(var(data_all(1:NoD*Ntri,:)-Ctot*S)/var(data_all(1:NoD*Ntri,:))))*100;
mean_tot_VAF_all = mean(tot_VAF);

%% NEW mean VAF based on Frere and Hug 2012, initially by Torres-Oviedo 2006
sum_sum_e2 = sum(sum((data_all-Ctot*S).^2));
sum_sum_E2 = sum(sum((data_all).^2));
tot_VAF_Torres = (1-sum_sum_e2/sum_sum_E2)*100;

%% Mean Synergy Contributions to Overall Reconstrcution
% for i = 1:NoS
%   tot_VAF_Torres_syn(i) =  (1-sum(sum((data_all-Ctot(:,i)*S(i,:)).^2))/sum_sum_E2)*100; 
%   tot_VAF_syn(i) =  (1-var((data_all-Ctot(:,i)*S(i,:)).^2)/var(data_all))*100;   
% tot_contr = sum(sum(contr(:,:,i)))/sum(sum(sum(contr(:,:,:))));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options, options2]=optim_Settings
options = optimset;
options = optimset(options,'Display', 'iter');
options = optimset(options,'Algorithm', 'sqp');
options = optimset(options,'TolFun', 1e-6);
options = optimset(options,'MaxFunEvals', 1e8);
options = optimset(options,'MaxIter', 1e8);
options = optimset(options,'TolX', 1e-15);
options = optimset(options,'TolCon', 1e-3);

options2 = optimset;
options2 = optimset(options2,'algorithm','active-set');
options2 = optimset(options2,'TolFun', 1e-6);
options2 = optimset(options2,'MaxFunEvals', 1e8);
options2 = optimset(options2,'MaxIter', 1e8);
options2 = optimset(options2,'TolX', 1e-15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function contr = PLOT_CNMF(wsize,Ntri,NoS,S,NoM,plot_sub_no,NoD,Ctot,data_all)

close all;
syn_col={'b','g','r','m','c','y'};
clear temp;
muscles={'GLTMAX','GLTMED','TFL','RFEM','VMED','VLAT','BFEM','STEN','TA','MGAS','LGAS','SOL'};
if isempty(findall(0,'type','figure'))
fh=0;
else
    fh=findall(0,'type','figure');
end
figure(fh+1);
set(figure(fh+1),'units','normalized','outerposition',[0 0 1*wsize 1]);
for j=1:Ntri
    if j==1

for i=1:NoS
        
temp=S;
    subplot(1,NoS,i)
    set(gca,'fontsize',24)

barh(temp(i,:),0.5,syn_col{i});hold all

ylim([0 NoM+1])

box on; grid on;
if i==1
set(gca, 'YTickLabel', muscles);
else
    set(gca, 'YTickLabel', []);
end
title(['S' num2str(i)])
end
   end
end

%% Coefficients 
fh=max(findall(0,'type','figure').Number);
for kkk=1:ceil(Ntri/plot_sub_no)
figure(kkk+fh)
set(figure(kkk+fh),'units','normalized','outerposition',[0 0 1*wsize 1]);

idx={'1','2','3','4','5'};
for j=(kkk-1)*plot_sub_no+1:min(kkk*plot_sub_no,Ntri) %j=1:Ntri
temp=Ctot((j-1)*NoD+1:j*NoD,:);
C_sub{j}=temp;
    for i=1:NoS
                
%         subplot(Ntri,NoS,i+NoS*(j-1))
subplot(plot_sub_no,NoS,i+(j-1-(kkk-1)*plot_sub_no)*NoS); 
        set(gca,'fontsize',24)

            plot(0:1:NoD-1,temp(:,i),syn_col{i},'LineWidth',3);

        set(gca,'Ylim',[0 100])
        set(gca,'Xlim',[0 NoD])
        box on; grid on;
        
        if j==1
            title(['C' num2str(i)]);
        end
        
        if j~=min(kkk*plot_sub_no,Ntri)
            set(gca,'XTickLabel',[])
        end
        if i==1
           set(gca, 'YTick', [0,100]); 
           ylabel(['trial', num2str(j)]);
        else
            set(gca,'YTickLabel',[]);
        end
        set(gca, 'XTick', 0:NoD/2:NoD);
    end
    
end
end

%% Reconstruction PLots
fh=max(findall(0,'type','figure').Number);
for kkk=1:ceil(Ntri/plot_sub_no)
    figure(kkk+fh)
set(figure(kkk+fh),'units','normalized','outerposition',[0 0 1 1]);
    
    for j=(kkk-1)*plot_sub_no+1:min(kkk*plot_sub_no,Ntri)
        A=data_all(1+NoD*(j-1):NoD*j,:);
    for i=1:NoM


    C=Ctot((j-1)*NoD+1:j*NoD,:);

subplot(plot_sub_no,NoM,i+(j-1-(kkk-1)*plot_sub_no)*NoM); 

    plot(0:NoD-1,A(:,i),'k',0:NoD-1,C*S(:,i),'--k','LineWidth',2);hold all;
    grid on;
    [r2(i+(j-1)*NoM),rmse(i+(j-1)*NoM)]=rsquare(A(:,i),C*S(:,i),false);

    
    ICC_recons(i+(j-1)*NoM)=ICC([A(:,i),C*S(:,i)], 'C-k', 0.05);
if max(A(:,i))<=0.1
   flag='NDom';
   r2_Dom(i+(j-1)*NoM)=NaN;
else
    flag='Dom';
    r2_Dom(i+(j-1)*NoM)=r2(i+(j-1)*NoM);
end
%             title(['r^2=',num2str(r2(i+(j-1)*NoM)),'   ','rmse=',num2str(rmse(i+(j-1)*NoM))],'Fontsize',12);
title({[flag,', r^2=',num2str(r2(i+(j-1)*NoM),2)];['ICC=',num2str(ICC_recons(i+(j-1)*NoM),2)]},'Fontsize',12);
            data_minmax(:,i+(j-1)*NoM)=[min(A(:,i));max(A(:,i))];
    for k=1:NoS
            plot(0:NoD-1,C(:,k)*S(k,i),syn_col{k},'LineWidth',2);
            if i==4
%                 legend('EMG','EMG_R','syn1','syn2','syn3','syn4','syn5')
                
            end
            %% contr of a size of Ntri (24 or 27) by NoM (12) by NoS (3, 4, or 5)
            contr(j,i,k)=trapz(1:NoD, C(:,k)*S(k,i))/trapz(1:NoD,C*S(:,i)); % e.g. the contribution of third synergy 
%             on muscle 6, trial 10 is contr(10,6,3).
%             Obviously, if there are 4 synergies, contr(10,6,1)+contr(10,6,2)+contr(10,6,3)+contr(10,6,4) =1
            
    end
    set(gca,'fontsize',12,'XTickLabel', [],'YTickLabel', [],'Ylim',[0,100],'Xlim',[0,NoD]);
    if i==1
       ylabel(['trial ', num2str(j)]); 
        set(gca,'YTick', [0,100],'YTickLabel', [0,100]);
    end
    if j==Ntri
       set(gca,'XTick', 0:NoD:NoD,'XTickLabel', 0:NoD:NoD);
    end
    if j==1
%        title({muscles{i};['r^2=',num2str(r2(i+(j-1)*NoM)),'   ','rmse=',num2str(rmse(i+(j-1)*8))]},'Fontsize',12) ;
title({muscles{i};[flag,', r^2=',num2str(r2(i+(j-1)*NoM),2)];['ICC=',num2str(ICC_recons(i+(j-1)*NoM),2)]},'Fontsize',12) ;
       
    end
    
    end
    end
    
end
    
%  
% figure(2);
% set(figure(2),'units','normalized','outerposition',[0 0 1*wsize 1]);
% 
% idx={'1','2','3','4','5'};
% for j=1:Ntri
% temp=Ctot((j-1)*NoD+1:j*NoD,:);
% C_sub{j}=temp;
%     for i=1:NoS
%                 
%         subplot(Ntri,NoS,i+NoS*(j-1))
%         set(gca,'fontsize',24)
% 
%             plot(0:1:NoD-1,temp(:,i),syn_col{i},'LineWidth',3);
% 
%         set(gca,'Ylim',[0 100])
%         set(gca,'Xlim',[0 NoD])
%         box on; grid on;
%         
%         if j==1
%             title(['C' num2str(i)]);
%         end
%         
%         if j~=Ntri
%             set(gca,'XTickLabel',[])
%         end
%         if i==1
%            set(gca, 'YTick', 0:100); 
%         else
%             set(gca,'YTickLabel',[]);
%         end
%         set(gca, 'XTick', 0:25:NoD);
%     end
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, table] = anova_rm(X, displayopt)
%   [p, table] = anova_rm(X, displayopt)
%   Single factor, repeated measures ANOVA.
%
%   [p, table] = anova_rm(X, displayopt) performs a repeated measures ANOVA
%   for comparing the means of two or more columns (time) in one or more
%   samples(groups). Unbalanced samples (i.e. different number of subjects 
%   per group) is supported though the number of columns (followups)should 
%   be the same. 
%
%   DISPLAYOPT can be 'on' (the default) to display the ANOVA table, or 
%   'off' to skip the display. For a design with only one group and two or 
%   more follow-ups, X should be a matrix with one row for each subject. 
%   In a design with multiple groups, X should be a cell array of matrixes.
% 
%   Example: Gait-Cycle-times of a group of 7 PD patients have been
%   measured 3 times, in one baseline and two follow-ups:
%
%   patients = [
%    1.1015    1.0675    1.1264
%    0.9850    1.0061    1.0230
%    1.2253    1.2021    1.1248
%    1.0231    1.0573    1.0529
%    1.0612    1.0055    1.0600
%    1.0389    1.0219    1.0793
%    1.0869    1.1619    1.0827 ];
%
%   more over, a group of 8 controls has been measured in the same protocol:
%
%   controls = [
%     0.9646    0.9821    0.9709
%     0.9768    0.9735    0.9576
%     1.0140    0.9689    0.9328
%     0.9391    0.9532    0.9237
%     1.0207    1.0306    0.9482
%     0.9684    0.9398    0.9501
%     1.0692    1.0601    1.0766
%     1.0187    1.0534    1.0802 ];
%
%   We are interested to see if the performance of the patients for the
%   followups were the same or not:
%  
%   p = anova_rm(patients);
%
%   By considering the both groups, we can also check to see if the 
%   follow-ups were significantly different and also check two see that the
%   two groups had a different performance:
%
%   p = anova_rm({patients controls});
%
%
%   ref: Statistical Methods for the Analysis of Repeated Measurements, 
%     C. S. Daivs, Springer, 2002
%
%   Copyright 2008, Arash Salarian
%   mailto://arash.salarian@ieee.org
%

if nargin < 2
    displayopt = 'on';
end

if ~iscell(X)
    X = {X};
end

%number of groups
s = size(X,2);  

%subjects per group 
n_h = zeros(s, 1);
for h=1:s
    n_h(h) = size(X{h}, 1);    
end
n = sum(n_h);

%number of follow-ups
t = size(X{1},2);   

% overall mean
y = 0;
for h=1:s
    y = y + sum(sum(X{h}));
end
y = y / (n * t);

% allocate means
y_h = zeros(s,1);
y_j = zeros(t,1);
y_hj = zeros(s,t);
y_hi = cell(s,1);
for h=1:s
    y_hi{h} = zeros(n_h(h),1);
end

% group means
for h=1:s
    y_h(h) = sum(sum(X{h})) / (n_h(h) * t);
end

% follow-up means
for j=1:t
    y_j(j) = 0;
    for h=1:s
        y_j(j) = y_j(j) + sum(X{h}(:,j));
    end
    y_j(j) = y_j(j) / n;
end

% group h and time j mean
for h=1:s
    for j=1:t
        y_hj(h,j) = sum(X{h}(:,j) / n_h(h));
    end
end

% subject i'th of group h mean
for h=1:s
    for i=1:n_h(h)
        y_hi{h}(i) = sum(X{h}(i,:)) / t;
    end
end

% calculate the sum of squares
ssG = 0;
ssSG = 0;
ssT = 0;
ssGT = 0;
ssR = 0;

for h=1:s
    for i=1:n_h(h)
        for j=1:t
            ssG  = ssG  + (y_h(h) - y)^2;
            ssSG = ssSG + (y_hi{h}(i) - y_h(h))^2;
            ssT  = ssT  + (y_j(j) - y)^2;
            ssGT = ssGT + (y_hj(h,j) - y_h(h) - y_j(j) + y)^2;
            ssR  = ssR  + (X{h}(i,j) - y_hj(h,j) - y_hi{h}(i) + y_h(h))^2;
        end
    end
end

% calculate means
if s > 1
    msG  = ssG  / (s-1);
    msGT = ssGT / ((s-1)*(t-1));
end
msSG = ssSG / (n-s);
msT  = ssT  / (t-1);
msR  = ssR  / ((n-s)*(t-1));


% calculate the F-statistics
if s > 1
    FG  = msG  / msSG;
    FGT = msGT / msR;
end
FT  = msT  / msR;
FSG = msSG / msR;


% single or multiple sample designs?
if s > 1
    % case for multiple samples
    pG  = 1 - fcdf(FG, s-1, n-s);
    pT  = 1 - fcdf(FT, t-1, (n-s)*(t-1));
    pGT = 1 - fcdf(FGT, (s-1)*(t-1), (n-s)*(t-1));
    pSG = 1 - fcdf(FSG, n-s, (n-s)*(t-1));

    p = [pT, pG, pSG, pGT];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Group' ssG s-1 msG FG pG
        'Ineratcion' ssGT (s-1)*(t-1) msGT FGT pGT
        'Subjects (matching)' ssSG n-s msSG FSG pSG
        'Error' ssR (n-s)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'multi-sample repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
else
    % case for only one sample
    pT  = 1 - fcdf(FT, t-1, (n-s)*(t-1));
    pSG = 1 - fcdf(FSG, n-s, (n-s)*(t-1));

    p = [pT, pSG];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Subjects (matching)' ssSG n-s msSG FSG pSG
        'Error' ssR (n-s)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
% Intraclass correlation
%   [r, LB, UB, F, df1, df2, p] = ICC(M, type, alpha, r0)
%
%   M is matrix of observations. Each row is an object of measurement and
%   each column is a judge or measurement.
%
%   'type' is a string that can be one of the six possible codes for the desired
%   type of ICC:
%       '1-1': The degree of absolute agreement among measurements made on
%         randomly seleted objects. It estimates the correlation of any two
%         measurements.
%       '1-k': The degree of absolute agreement of measurements that are
%         averages of k independent measurements on randomly selected
%         objects.
%       'C-1': case 2: The degree of consistency among measurements. Also known
%         as norm-referenced reliability and as Winer's adjustment for
%         anchor points. case 3: The degree of consistency among measurements maded under
%         the fixed levels of the column factor. This ICC estimates the
%         corrlation of any two measurements, but when interaction is
%         present, it underestimates reliability.
%       'C-k': case 2: The degree of consistency for measurements that are
%         averages of k independent measurements on randomly selected
%         onbjectgs. Known as Cronbach's alpha in psychometrics. case 3:  
%         The degree of consistency for averages of k independent
%         measures made under the fixed levels of column factor.
%       'A-1': case 2: The degree of absolute agreement among measurements. Also
%         known as criterion-referenced reliability. case 3: The absolute 
%         agreement of measurements made under the fixed levels of the column factor.
%       'A-k': case 2: The degree of absolute agreement for measurements that are
%         averages of k independent measurements on randomly selected objects.
%         case 3: he degree of absolute agreement for measurements that are
%         based on k independent measurements maded under the fixed levels
%         of the column factor.
%
%       ICC is the estimated intraclass correlation. LB and UB are upper
%       and lower bounds of the ICC with alpha level of significance. 
%
%       In addition to estimation of ICC, a hypothesis test is performed
%       with the null hypothesis that ICC = r0. The F value, degrees of
%       freedom and the corresponding p-value of the this test are
%       reported.
%
%       (c) Arash Salarian, 2008
%
%       Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About
%       Some Intraclass Correlation Coefficients", Psychological Methods,
%       Vol. 1, No. 1, pp. 30-46, 1996
%

if nargin < 3
    alpha = .05;
end

if nargin < 4
    r0 = 0;
end

[n, k] = size(M);
[p, table] = anova_rm(M, 'off');

SSR = table{3,2};
SSE = table{4,2};
SSC = table{2,2};
SSW = SSE + SSC;

MSR = SSR / (n-1);
MSE = SSE / ((n-1)*(k-1));
MSC = SSC / (k-1);
MSW = SSW / (n*(k-1));

switch type
    case '1-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case '1-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'C-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'C-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'A-1'
        [r, LB, UB, F, df1, df2, p] = ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k);
    case 'A-k'
        [r, LB, UB, F, df1, df2, p] = ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k);
end


%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSW) / (MSR + (k-1)*MSW);

F = (MSR/MSW) * (1-r0)/(1+(k-1)*r0);
df1 = n-1;
df2 = n*(k-1);
p = 1-fcdf(F, df1, df2);

FL = F / finv(1-alpha/2, n-1, n*(k-1));
FU = F * finv(1-alpha/2, n*(k-1), n-1);

LB = (FL - 1) / (FL + (k-1));
UB = (FU - 1) / (FU + (k-1));

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSW) / MSR;

F = (MSR/MSW) * (1-r0);
df1 = n-1;
df2 = n*(k-1);
p = 1-fcdf(F, df1, df2);

FL = F / finv(1-alpha/2, n-1, n*(k-1));
FU = F * finv(1-alpha/2, n*(k-1), n-1);

LB = 1 - 1 / FL;
UB = 1 - 1 / FU;

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (k-1)*MSE);

F = (MSR/MSE) * (1-r0)/(1+(k-1)*r0);
df1 = n - 1;
df2 = (n-1)*(k-1);
p = 1-fcdf(F, df1, df2);

FL = F / finv(1-alpha/2, n-1, (n-1)*(k-1));
FU = F * finv(1-alpha/2, (n-1)*(k-1), n-1);

LB = (FL - 1) / (FL + (k-1));
UB = (FU - 1) / (FU + (k-1));

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k) %% Double-Checked...Moh
r = (MSR - MSE) / MSR;

F = (MSR/MSE) * (1-r0);
df1 = n - 1;
df2 = (n-1)*(k-1); 
p = 1-fcdf(F, df1, df2);

FL = F / finv(1-alpha/2, n-1, (n-1)*(k-1));
FU = F * finv(1-alpha/2, (n-1)*(k-1), n-1);

LB = 1 - 1 / FL;
UB = 1 - 1 / FU;

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (k-1)*MSE + k*(MSC-MSE)/n);

a = (k*r0) / (n*(1-r0));
b = 1 + (k*r0*(n-1))/(n*(1-r0));
F = MSR / (a*MSC + b*MSE);
df1 = n - 1;
df2 = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));
p = 1-fcdf(F, df1, df2);

a = k*r/(n*(1-r));
b = 1+k*r*(n-1)/(n*(1-r));
v = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));

Fs = finv(1-alpha/2, n-1, v);
LB = n*(MSR - Fs*MSE)/(Fs*(k*MSC + (k*n - k - n)*MSE) + n*MSR);

Fs = finv(1-alpha/2, v, n-1);
UB = n*(Fs*MSR-MSE)/(k*MSC + (k*n - k - n)*MSE + n*Fs*MSR);

%----------------------------------------
function [r, LB, UB, F, df1, df2, p] = ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k)
r = (MSR - MSE) / (MSR + (MSC-MSE)/n);

c = r0/(n*(1-r0));
d = 1 + (r0*(n-1))/(n*(1-r0));
F = MSR / (c*MSC + d*MSE);
df1 = n - 1;
df2 = (c*MSC + d*MSE)^2/((c*MSC)^2/(k-1) + (d*MSE)^2/((n-1)*(k-1)));
p = 1-fcdf(F, df1, df2);

a = r/(n*(1-r));
b = 1+r*(n-1)/(n*(1-r));
v = (a*MSC + b*MSE)^2/((a*MSC)^2/(k-1) + (b*MSE)^2/((n-1)*(k-1)));

Fs = finv(1-alpha/2, n-1, v);
LB = n*(MSR - Fs*MSE)/(Fs*(MSC-MSE) + n*MSR);

Fs = finv(1-alpha/2, v, n-1);
UB = n*(Fs*MSR - MSE)/(MSC - MSE + n*Fs*MSR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r2 rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));
