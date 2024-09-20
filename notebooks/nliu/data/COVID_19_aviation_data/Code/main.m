

% Adjacency matrix, flight matrix A1, A2, A3


addpath  ../Pre_processing/COVID19_datasets

stateCode=[9 25 34 36 44 ];

counties_all=uscounties_confirmed.FIPS;

counties = grab_County_fips(stateCode,counties_all);

countyadjacency2010=readtable('county_adjacency2010.csv');
countyIndices=countyadjacency2010.fipscounty;
countyAdjacency=countyadjacency2010.fipsneighbor;

A = grab_Adjacency(stateCode,counties,countyIndices,countyAdjacency);

A2 = eye(length(counties));
A1 =A-A2;

addpath ../Data/Flight_outputs
load('A3_profile.mat')

current_data_date=datenum (2020,11,18);
flight_update_date= datenum(2020,7,31);

days=current_data_date-flight_update_date;

A3_last=repmat(A3(:,:,end),[1 1 days]); % repeat the last day of Jul to make up the matrix size

A3=cat(3, A3, A3_last);


csvwrite('../Data/SEIR_sim_outputs/county_adjacency.csv', A1);


%% p, r, e 
ind=find(ismember(counties_all,counties));
sel_cty=den_cat(ind,:);
% csvwrite('County_density_list.csv',sel_cty)

ind_ub_den=ind(find(sel_cty(:,3)==1));
ind_rr_den=ind(find(sel_cty(:,3)==2));

A_ub_ind=find(sel_cty(:,3)==1);
A_rr_ind=find(sel_cty(:,3)==2);

data_confirmed_county=norm_confirmed(:,ind);

% data_deaths_county=norm_deaths(:,ind);
% data_recovered_county=norm_recovered(:,ind);
% data_infected=data_confirmed_county-data_deaths_county-data_recovered_county;
% data_removed=data_deaths_county+data_recovered_county;
% 
% 
% for q=2:length(data_recovered_county(:,1))
%     for qt=1:length (data_recovered_county(1,:))
%         %data_recovered_county(q,qt)==0 | 
%     if data_recovered_county(q,qt)<data_recovered_county(q-1,qt)
%         data_recovered_county(q,qt)=data_recovered_county(q-1,qt);
%     end
%     end
% end

[roww,coll]=size(data_confirmed_county);
data_exposed =zeros (roww,coll);

%removed data
pre_removed=21;

for ii=(pre_removed+1):roww
data_removed (ii,:)=data_confirmed_county(ii-pre_removed,:);
end
data_removed(data_removed<0)=0; 

%infected data
data_infected =data_confirmed_county -data_removed ;

data_infected(data_infected <0)=0; 

% exposed data 


pre_exposed=14; % can be changed 

for ii=1:roww-pre_exposed
 data_exposed(ii,:)=data_confirmed_county(ii+pre_exposed,:);
end

data_exposed(data_exposed <0)=0; 

 
%% Simulation initial time

start_date= datetime(2020,3,20);
end_date= datetime(2020,8,4);

start_timebx= datevec(start_date);
start_timeok= start_timebx(1)*10000+ start_timebx(2)*100+ start_timebx(3);
[r5row, r5col]=find(recover5==start_timeok);
start_time=r5col-1;
simtime=365;
h = 1;



 


%% Special case (eliminate zero columns)
 
zero_col=find(all(data_confirmed_county==0));
zero_county=counties(zero_col);
data_infected(:, zero_col)=[];
data_removed(:, zero_col)=[];
data_exposed(:, zero_col)=[];

A1(:, zero_col)=[];
A1(zero_col, :)=[];

A2(:, zero_col)=[];
A2(zero_col, :)=[];

A3(:, zero_col,:)=[];
A3(zero_col, :,:)=[];

sel_cty(zero_col,:)=[];
A_ub_ind=find(sel_cty(:,3)==1);
A_rr_ind=find(sel_cty(:,3)==2);







%% SEIR regional


close all



[beta_ub,beta_rr,betae_ub,betae_rr, sigma,gamma, betaGammaMu] = Calibration_SEIR(data_infected(1:end-pre_exposed,:),data_removed(1:end-pre_exposed,:),data_exposed(1:end-pre_exposed,:),h,A1,A2,A3, pre_exposed,pre_removed,A_ub_ind,A_rr_ind); % learned
%  sigma=[0.3;0.1];
[sim_ub,rsim_ub,esim_ub, sim_rr,rsim_rr,esim_rr]= Simulation_SEIR(h,betae_ub, betae_rr, beta_ub, beta_rr,A1,A2,A3,data_infected(start_time,:),simtime,data_removed(start_time,:),data_exposed(start_time,:),gamma,sigma, A_ub_ind, A_rr_ind);

csvwrite('../Data/SEIR_sim_outputs/sim_infected_level_urban.csv', sim_ub);
csvwrite('../Data/SEIR_sim_outputs/sim_infected_level_rural.csv', sim_rr);
realstart=datetime(2020,1,22); % data start time 
real_date = realstart + caldays(0:length(data_infected(:,1))-1);
realend=real_date(end); % data end time 


box on


data_infected_ub= data_infected(:, A_ub_ind);
data_infected_rr= data_infected(:, A_rr_ind);

data_removed_ub= data_removed(:, A_ub_ind);
data_removed_rr= data_removed(:, A_rr_ind);

data_exposed_ub= data_exposed(:, A_ub_ind);
data_exposed_rr= data_exposed(:, A_rr_ind);

date = start_date + caldays(0:simtime-1);
bg=find(real_date==start_date);
ed=find(real_date==realend);
x_real_ub=data_infected_ub(bg:ed,:);
x_real_rr=data_infected_rr(bg:ed,:);
date_real = start_date + caldays(0:ed-bg);
x_sim_ub=sim_ub(1:(ed-bg+1),:);
x_sim_rr=sim_rr(1:(ed-bg+1),:);

r_real_ub=data_removed_ub(bg:ed,:);
r_real_rr=data_removed_rr(bg:ed,:);
r_sim_ub=rsim_ub(1:(ed-bg+1),:);
r_sim_rr=rsim_rr(1:(ed-bg+1),:);

e_real_ub=data_exposed_ub(bg:ed-pre_exposed,:);
e_real_rr=data_exposed_rr(bg:ed-pre_exposed,:);
e_sim_ub=esim_ub(1:(ed-bg+1-pre_exposed),:);
e_sim_rr=esim_rr(1:(ed-bg+1-pre_exposed),:);

% figure(1)
% plot(date_real,x_sim_ub,'LineWidth',2)
% 
% set(gca,'Fontsize',14)
% 
% ylabel('Infection level x(t)')
% 
% title('Simulated x(t): Urban')
% set(gca,'Fontsize',14)
% outFile = sprintf('Figure/SEIR_ub.png');
% saveas (gcf, outFile)

figure(1)

newcolors = [0 0.4470 0.7410
            0.8500 0.3250 0.0980
            0.9290 0.6940 0.1250
            0.4940 0.1840 0.5560
            0.4660 0.6740 0.1880
            0.3010 0.7450 0.9330
            0.6350 0.0780 0.1840];
            
newcolors=repmat(newcolors,10,1);

[days, num_rr]=size(x_sim_rr);
[days, num_ub]=size(x_sim_ub);

for i=1:num_rr

plot(date_real,x_sim_rr(:,i),'Color',newcolors(i,:),'LineWidth',2)

hold on
plot(date_real,x_real_rr(:,i),'-.','Color',newcolors(i,:),'LineWidth',1)

%ylabel('p(t)')
%title('Simulated x(t): Rural')

xlim([start_date end_date]);
ylim([0 0.025])
pbaspect([3 1 1])
set(gca,'Fontsize',18)
end

outFile = sprintf('../Data/SEIR_sim_outputs/Figure/SEIR_rr_w.png');
saveas (gcf, outFile)

[days, num_ub]=size(x_sim_ub);

figure(2)

for i=1:num_ub
plot(date_real,x_sim_ub(:,i),'Color',newcolors(i,:),'LineWidth',2)
hold on
plot(date_real,x_real_ub(:,i),'-.','Color',newcolors(i,:),'LineWidth',1)

end

%ylabel('p(t)')
%title('Simulated x(t): Urban')
xlim([start_date end_date]);
ylim([0 0.025])
pbaspect([3 1 1])
set(gca,'Fontsize',18)
addpath ../Data/SEIR_sim_outputs
outFile = sprintf('../Data/SEIR_sim_outputs/Figure/SEIR_ub_w.png');
saveas (gcf, outFile)

nm=2;
x_error_ub = x_real_ub-x_sim_ub;
xn_err_ub= (norm(x_error_ub,nm))/(norm(x_real_ub,nm))

x_error_rr = x_real_rr-x_sim_rr;
xn_err_rr= (norm(x_error_rr,nm))/(norm(x_real_rr,nm))

nm=2;
r_error_ub = r_real_ub-r_sim_ub;
rn_err_ub= (norm(r_error_ub,nm))/(norm(r_real_ub,nm))

r_error_rr = r_real_rr-r_sim_rr;
rn_err_rr= (norm(r_error_rr,nm))/(norm(r_real_rr,nm))


nm=2;
e_error_ub = e_real_ub-e_sim_ub;
en_err_ub= (norm(e_error_ub,nm))/(norm(e_real_ub,nm))

e_error_rr = e_real_rr-e_sim_rr;
en_err_rr= (norm(e_error_rr,nm))/(norm(e_real_rr,nm))











