
close all
addpath ../Pre_processing/Transportation_air_network_datasets
Jul20= readtable("Jul200.csv");
Jun20= readtable("Jun200.csv");
May20= readtable("May200.csv");
Apr20= readtable("Apr200.csv");
Mar20= readtable("Mar200.csv");
Feb20= readtable("Feb200.csv");
Jan20= readtable("Jan200.csv");
fleet= readtable("fleet.csv");
% addpath ../Pre_processing/Transportation_air_network_datasets
tsf=readtable("../Pre_processing/Transportation_air_network_datasets/Airport_fips_northeastern.csv");

Airlines={'AA';'UA';'DL';'AS';'NK';'B6';'F9';'WN';'MQ';'OH';'OO';'YV';'YX';'9E';'EV'};% Select your airlines
State={'Connecticut','Massachusetts','New Jersy','New York','Rhode Island'}; %select your state 

stateCode=[9 25 34 36 44];

counties_all=uscounties_confirmed.FIPS; % Grab the counties in these five states, run 'COVID_processing' first to access uscounties_confirmed.FIPS 

addpath ../Code
indices = grab_County_fips(stateCode,counties_all); % fips of selected counties
%%
% Jul
cJul20=table2cell(Jul20);
[Lia, Lcob]=ismember(cJul20(:,4),Airlines);
first_usefulJul20=cJul20(Lia,:);
[Lia, Lcob]=ismember(first_usefulJul20(:,9),State);
second_usefulJul20=first_usefulJul20(Lia,:);
[Lia, Lcob]=ismember(second_usefulJul20(:,12),State);
third_usefulJul20=second_usefulJul20(Lia,:);

[Lia, Lcob]=ismember(third_usefulJul20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJul20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJul20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJul20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJul20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulJul20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJul20(i,16)=third_usefulJul20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulJul20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJul20(i,16)=third_usefulJul20(i+1,16);
    end
end


t=31;

nn=length(indices);

A3_Jul=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulJul20(:,2))==i);
    flight_day=third_usefulJul20(tloc,:);
    dep_arr_Jul=flight_day(:,14:16);
    mat_dep_arr_Jul=cell2mat(dep_arr_Jul);
    [C,ia,ic]=unique(mat_dep_arr_Jul(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Jul(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Jul(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
    
end

% Jun
cJun20=table2cell(Jun20);
[Lia, Lcob]=ismember(cJun20(:,4),Airlines);
first_usefulJun20=cJun20(Lia,:);
[Lia, Lcob]=ismember(first_usefulJun20(:,9),State);
second_usefulJun20=first_usefulJun20(Lia,:);
[Lia, Lcob]=ismember(second_usefulJun20(:,12),State);
third_usefulJun20=second_usefulJun20(Lia,:);

[Lia, Lcob]=ismember(third_usefulJun20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJun20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJun20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJun20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJun20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulJun20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJun20(i,16)=third_usefulJun20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulJun20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJun20(i,16)=third_usefulJun20(i+1,16);
    end
end


t=30;

nn=length(indices);

A3_Jun=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulJun20(:,2))==i);
    flight_day=third_usefulJun20(tloc,:);
    dep_arr_Jun=flight_day(:,14:16);
    mat_dep_arr_Jun=cell2mat(dep_arr_Jun);
    [C,ia,ic]=unique(mat_dep_arr_Jun(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Jun(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Jun(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
    
end

% May 
cMay20=table2cell(May20);
[Lia, Lcob]=ismember(cMay20(:,4),Airlines);
first_usefulMay20=cMay20(Lia,:);
[Lia, Lcob]=ismember(first_usefulMay20(:,9),State);
second_usefulMay20=first_usefulMay20(Lia,:);
[Lia, Lcob]=ismember(second_usefulMay20(:,12),State);
third_usefulMay20=second_usefulMay20(Lia,:);

[Lia, Lcob]=ismember(third_usefulMay20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulMay20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulMay20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulMay20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulMay20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulMay20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulMay20(i,16)=third_usefulMay20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulMay20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulMay20(i,16)=third_usefulMay20(i+1,16);
    end
end


t=31;

nn=length(indices);

A3_May=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulMay20(:,2))==i);
    flight_day=third_usefulMay20(tloc,:);
    dep_arr_May=flight_day(:,14:16);
    mat_dep_arr_May=cell2mat(dep_arr_May);
    [C,ia,ic]=unique(mat_dep_arr_May(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_May(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_May(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
    
end

% Apr 
cApr20=table2cell(Apr20);
[Lia, Lcob]=ismember(cApr20(:,4),Airlines);
first_usefulApr20=cApr20(Lia,:);
[Lia, Lcob]=ismember(first_usefulApr20(:,9),State);
second_usefulApr20=first_usefulApr20(Lia,:);
[Lia, Lcob]=ismember(second_usefulApr20(:,12),State);
third_usefulApr20=second_usefulApr20(Lia,:);

[Lia, Lcob]=ismember(third_usefulApr20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulApr20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulApr20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulApr20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulApr20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulApr20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulApr20(i,16)=third_usefulApr20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulApr20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulApr20(i,16)=third_usefulApr20(i+1,16);
    end
end


t=30;

nn=length(indices);

A3_Apr=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulApr20(:,2))==i);
    flight_day=third_usefulApr20(tloc,:);
    dep_arr_Apr=flight_day(:,14:16);
   mat_dep_arr_Apr=cell2mat(dep_arr_Apr);
    [C,ia,ic]=unique(mat_dep_arr_Apr(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Apr(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Apr(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
    
end

% Mar 
cMar20=table2cell(Mar20);
[Lia, Lcob]=ismember(cMar20(:,4),Airlines);
first_usefulMar20=cMar20(Lia,:);
[Lia, Lcob]=ismember(first_usefulMar20(:,9),State);
second_usefulMar20=first_usefulMar20(Lia,:);
[Lia, Lcob]=ismember(second_usefulMar20(:,12),State);
third_usefulMar20=second_usefulMar20(Lia,:);

[Lia, Lcob]=ismember(third_usefulMar20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulMar20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulMar20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulMar20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulMar20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulMar20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulMar20(i,16)=third_usefulMar20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulMar20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulMar20(i,16)=third_usefulMar20(i+1,16);
    end
end


t=31;

nn=length(indices);

A3_Mar=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulMar20(:,2))==i);
    flight_day=third_usefulMar20(tloc,:);
    dep_arr_Mar=flight_day(:,14:16);
    mat_dep_arr_Mar=cell2mat(dep_arr_Mar);
    [C,ia,ic]=unique(mat_dep_arr_Mar(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Mar(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Mar(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
end

% Feb
cFeb20=table2cell(Feb20);
[Lia, Lcob]=ismember(cFeb20(:,4),Airlines);
first_usefulFeb20=cFeb20(Lia,:);
[Lia, Lcob]=ismember(first_usefulFeb20(:,9),State);
second_usefulFeb20=first_usefulFeb20(Lia,:);
[Lia, Lcob]=ismember(second_usefulFeb20(:,12),State);
third_usefulFeb20=second_usefulFeb20(Lia,:);

[Lia, Lcob]=ismember(third_usefulFeb20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulFeb20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulFeb20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulFeb20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulFeb20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulFeb20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulFeb20(i,16)=third_usefulFeb20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulFeb20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulFeb20(i,16)=third_usefulFeb20(i+1,16);
    end
end


t=29;

nn=length(indices);

A3_Feb=zeros(nn,nn,t);

for i=1:t
    tloc=find(cell2mat(third_usefulFeb20(:,2))==i);
    flight_day=third_usefulFeb20(tloc,:);
    dep_arr_Feb=flight_day(:,14:16);
    mat_dep_arr_Feb=cell2mat(dep_arr_Feb);
    [C,ia,ic]=unique(mat_dep_arr_Feb(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Feb(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Feb(Lorg,Ldes,i)=sum_dep_arr(k,4);
    end
    
    
    
end

% Jan
cJan20=table2cell(Jan20);
[Lia, Lcob]=ismember(cJan20(:,4),Airlines);
first_usefulJan20=cJan20(Lia,:);
[Lia, Lcob]=ismember(first_usefulJan20(:,9),State);
second_usefulJan20=first_usefulJan20(Lia,:);
[Lia, Lcob]=ismember(second_usefulJan20(:,12),State);
third_usefulJan20=second_usefulJan20(Lia,:);

[Lia, Lcob]=ismember(third_usefulJan20(:,7),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJan20(i,14)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJan20(:,10),table2cell(tsf(1:30,10)));
for i=1:length(Lcob)
third_usefulJan20(i,15)=table2cell(tsf(Lcob(i),11));
end

[Lia, Lcob]=ismember(third_usefulJan20(:,5),table2cell(fleet(:,7)));

for i=2:length(Lcob)
    if Lcob(i)~=0
        third_usefulJan20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJan20(i,16)=third_usefulJan20(i-1,16);
    end
end

for i=1
    if Lcob(i)~=0
        third_usefulJan20(i,16)=table2cell(fleet(Lcob(i),10));
    else 
        third_usefulJan20(i,16)=third_usefulJan20(i+1,16);
    end
end


t=31;

nn=length(indices);

A3_Jan=zeros(nn,nn,t-21);

for i=22:t
    tloc=find(cell2mat(third_usefulJan20(:,2))==i);
    flight_day=third_usefulJan20(tloc,:);
    dep_arr_Jan=flight_day(:,14:16);
    mat_dep_arr_Jan=cell2mat(dep_arr_Jan);
    [C,ia,ic]=unique(mat_dep_arr_Jan(:,1:2),'rows');
    a_sum = accumarray(ic, mat_dep_arr_Jan(:,3));
    sum_dep_arr=[C a_sum];
    sum_dep_arr(:,4)=a_sum/max(a_sum);
    for k=1:length(sum_dep_arr(:,1))
    Lorg=find(sum_dep_arr(k,1)==indices);
    Ldes=find(sum_dep_arr(k,2)==indices);
    A3_Jan(Lorg,Ldes,i-21)=sum_dep_arr(k,4);
    end
    
    
end

%% Generate A3 (flight matrix for every day from Jan-Jul)
% A3_end=repmat(A3_Jul(:,:,31),[1 1 4]);
A3=cat(3, A3_Jan, A3_Feb, A3_Mar, A3_Apr, A3_May,A3_Jun, A3_Jul);
save('../Data/A3_profile','A3') % A3 is a 110*110*161(time ) matrix 
csvwrite('../Data/Air_network_adjacency.csv', A3);

%% Generate flight info list
old_title=Jun20.Properties.VariableNames;
title=cat(2, old_title, { 'ORIGIN_FIPS','DES_FIPS','NUM_SEATS' });
third_useful= cell2table([third_usefulJan20;third_usefulFeb20;third_usefulMar20;third_usefulApr20;third_usefulMay20;third_usefulJun20;third_usefulJul20]);
third_useful.Properties.VariableNames=title;
writetable(third_useful,'../Data/List_of_flight.csv') 

