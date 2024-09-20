function [beta_ub,beta_rr,betae_ub,betae_rr, sigma,gamma, betaGammaMu] = mfID_SEIRA3_regional(sol,re,e,h,A1,A2,A3,pre1,pre,A_ub_ind,A_rr_ind)

%   sol=data;
%   re=data_r;
% e=data_e;
% together 
[m, n] = size(sol);

x = reshape(sol',m*n,1);
y = reshape(re',m*n,1);
z = reshape(e',m*n,1);
 
n_ub=length(A_ub_ind);
n_rr=length(A_rr_ind);

ub_sol=reshape(sol(:,A_ub_ind)',m*n_ub,1);
rr_sol=reshape(sol(:,A_rr_ind)',m*n_rr,1);
ub_re=reshape(re(:,A_ub_ind)',m*n_ub,1);
rr_re=reshape(re(:,A_rr_ind)',m*n_rr,1);
ub_e=reshape(e(:,A_ub_ind)',m*n_ub,1);
rr_e=reshape(e(:,A_rr_ind)',m*n_rr,1);

q=8;

Alpha = zeros(m*n-n,q);
phi= zeros(m*n-n,q);
Sigma = zeros(m*n-n,q);

Alpha_ub = zeros(m*n_ub-n_ub,q);
phi_ub= zeros(m*n_ub-n_ub,q);
Sigma_ub= zeros(m*n_ub-n_ub,q);

Alpha_rr = zeros(m*n_rr-n_rr,q);
phi_rr= zeros(m*n_rr-n_rr,q);
Sigma_rr= zeros(m*n_rr-n_rr,q);



% weight=ones((m-1)*n*2,1);
% weight(29*n:30*n-1,1)=20*ones(n,1);
% weight(n*m+29*n:n*m+30*n-1,1)=20*ones(n,1);

% Figure out different regions



diff = x(n+1:end) - x(1:m*n-n);
ydiff= y(n+1:end) - y(1:m*n-n);
zdiff =z(n+1:end) - z(1:m*n-n);

ub_diff=ub_sol(length(A_ub_ind)+1:end) - ub_sol(1:m*length(A_ub_ind)-length(A_ub_ind));
rr_diff=rr_sol(length(A_rr_ind)+1:end) - rr_sol(1:m*length(A_rr_ind)-length(A_rr_ind));

ub_ydiff=ub_re(length(A_ub_ind)+1:end) - ub_re(1:m*length(A_ub_ind)-length(A_ub_ind));
rr_ydiff=rr_re(length(A_rr_ind)+1:end) - rr_re(1:m*length(A_rr_ind)-length(A_rr_ind));

ub_zdiff=ub_e(length(A_ub_ind)+1:end) - ub_e(1:m*length(A_ub_ind)-length(A_ub_ind));
rr_zdiff=rr_e(length(A_rr_ind)+1:end) - rr_e(1:m*length(A_rr_ind)-length(A_rr_ind));

pdiff=[zdiff;diff;ydiff];
ub_pdiff=[ub_zdiff;ub_diff;ub_ydiff];
rr_pdiff=[rr_zdiff;rr_diff;rr_ydiff];


A1_ub=A1(A_ub_ind,:);
A1_rr=A1(A_rr_ind,:);
A2_ub=A2(A_ub_ind,:);
A2_rr=A2(A_rr_ind,:);
A3_ub3=A3(A_ub_ind,:,:);
A3_rr3=A3(A_rr_ind,:,:);

rf1a1=1; % NY
rf1a2=1; % NY
rf1a3=1; % NY

% rf2a1=0.7; % MN
% rf2a2=0.3; % MN

% rf2a1=0.7; % tx
% rf2a2=0.2; % tx

%  rf2a1=0.4;  % MI
%  rf2a2=0.3;  % MI

% rf2a1=0.55; % NY
% rf2a2=0.45; % NY
% % ub 
%  rf2a1=0.45; % NY
%  rf2a2=0.45; % NY

% % rr
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

for i = 1:m-1
    xi = x(n*(i-1)+1:n*i);
    yi = y(n*(i-1)+1:n*i);
    zi = z(n*(i-1)+1:n*i);
    xi_ub = ub_sol(n_ub*(i-1)+1:n_ub*i);
    yi_ub = ub_re(n_ub*(i-1)+1:n_ub*i);
    zi_ub = ub_e(n_ub*(i-1)+1:n_ub*i);
    
    xi_rr = rr_sol(n_rr*(i-1)+1:n_rr*i);
    yi_rr = rr_re(n_rr*(i-1)+1:n_rr*i);
    zi_rr = rr_e(n_rr*(i-1)+1:n_rr*i);

if i >= 65 && i<=90 %139
      A1_ub=rf2a1*A1_ub;
      A2_ub=rf2a2*A2_ub;
      A3_ub=rf2a3*A3_ub3(:,:,i);
      
%       A1_rr=rf2a1*A1_rr;
%       A2_rr=rf2a2*A2_rr;
%       A3_rr=rf2a3*A3_rr;
      A1_rr=1*A1_rr;
      A2_rr=1*A2_rr;
      A3_rr=rf2a3*A3_rr3(:,:,i);
elseif i > 90 && i<=139
      A1_ub=rf2a1*A1_ub;
      A2_ub=rf2a2*A2_ub;
      A3_ub=rf2a3*A3_ub3(:,:,i);

       A1_rr=0.8*A1_rr;
       A2_rr=0.85*A2_rr;
       A3_rr=rf2a3*A3_rr3(:,:,i);
% elseif i >= 72
%       A1=rf2a1*A1;
%       A2=rf2a2*A2;
 else  
      A1_ub=rf1a1*A1_ub;
      A2_ub=rf1a2*A2_ub;
      A3_ub=rf1a3*A3_ub3(:,:,i);
      
      A1_rr=rf1a1*A1_rr;
      A2_rr=rf1a2*A2_rr; 
      A3_rr=rf1a3*A3_rr3(:,:,i);
end 
    
    Alpha_ub((n_ub*(i-1)+1):(n_ub*i),:) = [h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A1_ub*zi, h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A2_ub*zi, h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A3_ub*zi, h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A1_ub*xi, h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A2_ub*xi,h*(eye(n_ub) - diag(xi_ub) - diag(yi_ub)-diag(zi_ub))*A3_ub*xi ,-h*zi_ub, zeros(length(xi_ub),1)];
    Sigma_ub((n_ub*(i-1)+1):(n_ub*i),:) = [zeros(length(xi_ub),1),zeros(length(xi_ub),1), zeros(length(xi_ub),1),zeros(length(xi_ub),1),zeros(length(xi_ub),1),zeros(length(xi_ub),1), h*zi_ub, -h*xi_ub];
    phi_ub((n_ub*(i-1)+1):(n_ub*i),:) = [zeros(length(xi_ub),1),zeros(length(xi_ub),1),zeros(length(xi_ub),1),zeros(length(xi_ub),1),zeros(length(xi_ub),1), zeros(length(xi_ub),1),zeros(length(xi_ub),1),h*xi_ub];
    
    Alpha_rr((n_rr*(i-1)+1):(n_rr*i),:) = [h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A1_rr*zi, h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A2_rr*zi, h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A3_rr*zi, h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A1_rr*xi, h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A2_rr*xi ,  h*(eye(n_rr) - diag(xi_rr) - diag(yi_rr)-diag(zi_rr))*A3_rr*xi, -h*zi_rr, zeros(length(xi_rr),1)];
    Sigma_rr((n_rr*(i-1)+1):(n_rr*i),:) = [zeros(length(xi_rr),1),zeros(length(xi_rr),1), zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1), h*zi_rr, -h*xi_rr];
    phi_rr((n_rr*(i-1)+1):(n_rr*i),:) = [zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1),zeros(length(xi_rr),1), zeros(length(xi_rr),1), h*xi_rr];
    
end



P=[Alpha_ub, zeros(m*n_ub-n_ub,q); Sigma_ub, zeros(m*n_ub-n_ub,q); phi_ub, zeros(m*n_ub-n_ub,q); zeros(m*n_rr-n_rr,q),Alpha_rr;zeros(m*n_rr-n_rr,q),Sigma_rr; zeros(m*n_rr-n_rr,q),phi_rr];
% size(P)
% size(weight)


%betaGammaMu = pinv(P)*[ub_pdiff;rr_pdiff;rr_pdiff]
betaGammaMu = pinv(P)*[ub_pdiff;rr_pdiff]
%betaGammaMu = pinv(P)*pdiff;



cvx_begin
cvx_solver mosek

variables be1_ub be2_ub be3_ub be1_rr be2_rr be3_rr b1_ub b2_ub b3_ub b1_rr b2_rr b3_rr gam_ub gam_rr %sig_ub sig_rr %b1_rr b2_rr 
%variables  b1 b2 gam 
% b0>=0


b1_ub>=0
b2_ub>=0
b3_ub>=0

b1_rr>=0
b2_rr>=0
b3_rr>=0

be1_ub>=0
be2_ub>=0
be3_ub>=0

be1_rr>=0
be2_rr>=0
be3_rr>=0

% b1_rr>=0
% b2_rr>=0

 gam_ub>=0
% 
gam_ub<=1/h

 gam_rr>=0
% 
gam_rr<=1/h

% sig_ub>=0
% % 
% sig_ub<=1/h
% 
% sig_rr>=0
% % 
% sig_rr<=1/h

for j=1:length(A1_ub(:,1))
    sum(A1_ub(j,:))*b1_ub+sum(A2_ub(j,:))*b2_ub+sum(A3_ub(j,:))*b3_ub+sum(A1_ub(j,:))*be1_ub+sum(A2_ub(j,:))*be2_ub+sum(A3_ub(j,:))*be3_ub  <= 1/h
end


for jj=1:length(A1_rr(:,1))
    sum(A1_rr(jj,:))*b1_rr+sum(A2_rr(jj,:))*b2_rr+sum(A3_rr(jj,:))*b3_rr+sum(A1_rr(jj,:))*be1_rr+sum(A2_rr(jj,:))*be2_rr+sum(A3_rr(jj,:))*be3_rr <= 1/h
end

 sig_ub=1/pre1;
 sig_rr=1/pre1;

% gam_ub=1/pre;
% gam_rr=1/pre;
%minimize(norm((P(:,1)*b1_ub + P(:,2)*b2_ub+ P(:,3)*gam+P(:,4)*b1_rr+ P(:,5)*b2_rr+ P(:,6)*gam+P(:,7)*b1_rr+ P(:,8)*b2_rr+ P(:,9)*gam)-[ub_pdiff;rr_pdiff;rr_pdiff],2))
%  minimize(norm(weight.*((P(:,1)*b1 + P(:,2)*b2+ P(:,3)*gam)-pdiff),2))
%minimize(norm((P(:,1)*b1_ub + P(:,2)*b2_ub+ P(:,3)*gam+P(:,4)*b1_rr+ P(:,5)*b2_rr+ P(:,6)*gam)-[ub_pdiff;rr_pdiff],2))
minimize(norm(P*[be1_ub;be2_ub;be3_ub;b1_ub;b2_ub; b3_ub;sig_ub;gam_ub; be1_rr;be2_rr;be3_rr; b1_rr;b2_rr;b3_rr; sig_rr;gam_rr]-[ub_pdiff;rr_pdiff],2))
cvx_end

betae_ub=[be1_ub;be2_ub;be3_ub]
beta_ub=[b1_ub;b2_ub;b3_ub]

betae_rr=[be1_rr;be2_rr;be3_rr]
beta_rr=[b1_rr;b2_rr;b3_rr]

sigma=[sig_ub;sig_rr]
gamma=[gam_ub;gam_rr]
betaGammaMu
%gammaMu = pinv(phi)*ydiff;