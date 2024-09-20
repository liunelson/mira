% This is the discrete euler's method approx. of the mean field approx of 
% the continuous 2^n mode; 
function [sol_ub,rsol_ub,esol_ub, sol_rr,rsol_rr,esol_rr]= dis_SEIRA3_regional(h,betae_ub, betae_rr, beta_ub, beta_rr,A1,A2,A3,x0,tspan,r0,e0,gamma,sigma, A_ub_ind, A_rr_ind)

n = length(x0);
n_ub =length(A_ub_ind);
n_rr =length(A_rr_ind);

A3_rep=repmat(A3(:,:,end),[1 1 tspan]);
A3=cat(3, A3,A3_rep);

A1_ub=A1(A_ub_ind,:);
A1_rr=A1(A_rr_ind,:);
A2_ub=A2(A_ub_ind,:);
A2_rr=A2(A_rr_ind,:);
A3_ub3=A3(A_ub_ind,:,:);
A3_rr3=A3(A_rr_ind,:,:);



sol = zeros(tspan,n);
rsol = zeros(tspan,n);
esol = zeros(tspan,n);


sol_ub = zeros(tspan,n_ub);
rsol_ub = zeros(tspan,n_ub);
esol_ub=zeros(tspan,n_ub);


sol_rr = zeros(tspan,n_rr);
rsol_rr = zeros(tspan,n_rr);
esol_rr=zeros(tspan,n_rr);

% sol(1,:) = x0';
% rsol(1,:) = r0';
x0_ub=x0(A_ub_ind);
x0_rr=x0(A_rr_ind);


r0_ub=r0(A_ub_ind);
r0_rr=r0(A_rr_ind);

e0_ub=e0(A_ub_ind);
e0_rr=e0(A_rr_ind);

sol(1,:) = x0;
rsol(1,:) = r0;
esol(1,:) = e0;


sol_ub(1,:) = x0_ub;
rsol_ub(1,:) = r0_ub;
esol_ub(1,:) = e0_ub;


sol_rr(1,:) = x0_rr;
rsol_rr(1,:) = r0_rr;
esol_rr(1,:) = e0_rr;

rf1a1=1; % NY
rf1a2=1; % NY
rf1a3=1; % NY
% rf2a1=0.7; % MN
% rf2a2=0.3; % MN

% rf2a1=0.7; % tx
% rf2a2=0.2; % tx

%  rf2a1=0.4;  % MI
%  rf2a2=0.3;  % MI
% 
% rf2a1=0.55; % NY
% rf2a2=0.45; % NY

% % ua 
%  rf2a1=0.45; % NY
%  rf2a2=0.45; % NY
 
%  % uc
%  rf2a1=0.54; % NY
%  rf2a2=0.45; % NY
 
 % rr
 rf2a1=0.6; % NY
 rf2a2=0.8; % NY
 rf2a3=0.6; % NY

% rf2a1=1; 
% rf2a2=1;


% rf3a1=0.9;
% rf3a2=0.3;

for i = 1:tspan-1
    if i >=7 && i <= 32 %81
        A1_ub=rf2a1*A1_ub;
        A2_ub=rf2a2*A2_ub;
        A3_ub=rf2a3*A3_ub3(:,:,i+59);
      
%       A1_rr=rf2a1*A1_rr;
%       A2_rr=rf2a2*A2_rr;
%       A3_rr=rf2a3*A3_rr;
      A1_rr=1*A1_rr;
      A2_rr=1*A2_rr;
      A3_rr=rf2a3*A3_rr3(:,:,i+59);
 elseif i > 32 && i <=81
      A1_ub=rf2a1*A1_ub;
      A2_ub=rf2a2*A2_ub;
      A3_ub=rf2a3*A3_ub3(:,:,i+59);

       A1_rr=0.8*A1_rr;
       A2_rr=0.85*A2_rr;
       A3_rr=rf2a3*A3_rr3(:,:,i+59);
 else  
      A1_ub=rf1a1*A1_ub;
      A2_ub=rf1a2*A2_ub;
      A3_ub=rf1a3*A3_ub3(:,:,i+59);
      
      A1_rr=rf1a1*A1_rr;
      A2_rr=rf1a2*A2_rr; 
      A3_rr=rf1a3*A3_rr3(:,:,i+59);
    end
    X = diag(sol(i,:));
    R = diag(rsol(i,:));
    E = diag(esol(i,:));
    
    
    X_ub = diag(sol_ub(i,:));
    R_ub = diag(rsol_ub(i,:));
    E_ub = diag(esol_ub(i,:));
     
    X_rr = diag(sol_rr(i,:));
    R_rr = diag(rsol_rr(i,:));
    E_rr = diag(esol_rr(i,:));
    
%      sol(i+1,:) = (sol(i,:)' + h*( (eye(n)-X-R)*beta(1)*A1*(sol(i,:)')+(eye(n)-X-R)*beta(2)*A2*(sol(i,:)')-gamma*(sol(i,:)')))';
% %      sol(i+1,:) = (sol(i,:)' + h*( (eye(n)-X-R)*beta(1)*A1*(sol(i,:)')+(eye(n)-X-R)*beta(2)*A2*(sol(i,:)')-gamma*X))';
%     rsol(i+1,:)= (rsol(i,:)'+h*gamma*(sol(i,:)'))';

    esol_ub(i+1,:) = (esol_ub(i,:)' + h*( (eye(n_ub)-X_ub-R_ub-E_ub)*beta_ub(1)*A1_ub*(sol(i,:)')+(eye(n_ub)-X_ub-R_ub-E_ub)*beta_ub(2)*A2_ub*(sol(i,:)')+(eye(n_ub)-X_ub-R_ub-E_ub)*beta_ub(3)*A3_ub*(sol(i,:)')+(eye(n_ub)-X_ub-R_ub-E_ub)*betae_ub(1)*A1_ub*(esol(i,:)')+(eye(n_ub)-X_ub-R_ub-E_ub)*betae_ub(2)*A2_ub*(esol(i,:)') +(eye(n_ub)-X_ub-R_ub-E_ub)*betae_ub(3)*A3_ub*(esol(i,:)')- sigma(1)*(esol_ub(i,:)')))';   
    sol_ub(i+1,:) = (sol_ub(i,:)' + h*(sigma(1)* esol_ub(i,:)' -gamma(1)*sol_ub(i,:)'))';
    rsol_ub(i+1,:)= (rsol_ub(i,:)'+h*(gamma(1)*(sol_ub(i,:)')))';
    
    esol_rr(i+1,:) = (esol_rr(i,:)' + h*( (eye(n_rr)-X_rr-R_rr-E_rr)*beta_rr(1)*A1_rr*(sol(i,:)')+(eye(n_rr)-X_rr-R_rr-E_rr)*beta_rr(2)*A2_rr*(sol(i,:)')+(eye(n_rr)-X_rr-R_rr-E_rr)*beta_rr(3)*A3_rr*(sol(i,:)')+(eye(n_rr)-X_rr-R_rr-E_rr)*betae_rr(1)*A1_rr*(esol(i,:)')+(eye(n_rr)-X_rr-R_rr-E_rr)*betae_rr(2)*A2_rr*(esol(i,:)')+(eye(n_rr)-X_rr-R_rr-E_rr)*betae_rr(3)*A3_rr*(esol(i,:)') - sigma(2)*(esol_rr(i,:)')))';   
    sol_rr(i+1,:) = (sol_rr(i,:)' + h*(sigma(2)* esol_rr(i,:)' -gamma(2)*sol_rr(i,:)'))';
    rsol_rr(i+1,:)= (rsol_rr(i,:)'+h*(gamma(2)*(sol_rr(i,:)')))';
    
     for k=1:length(A_ub_ind)
        sol(i+1,A_ub_ind(k))= sol_ub(i+1,k);
        esol(i+1,A_ub_ind(k))= esol_ub(i+1,k);
        rsol(i+1,A_ub_ind(k))= rsol_ub(i+1,k);
    end
    
     for k=1:length(A_rr_ind)
        sol(i+1,A_rr_ind(k))= sol_rr(i+1,k);
        esol(i+1,A_rr_ind(k))= esol_rr(i+1,k);
        rsol(i+1,A_rr_ind(k))= rsol_rr(i+1,k);
    end
 
end
