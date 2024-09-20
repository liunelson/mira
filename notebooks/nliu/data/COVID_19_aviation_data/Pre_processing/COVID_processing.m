clear all

close all
addpath ../Pre_processing/COVID19_datasets
% If you do not have source data file, please go to "Source file" folder
% for downloading. 
confirmedUS=readtable("time_series_covid19_confirmed_US.csv",'ReadVariableNames',1);
deathsUS=readtable("time_series_covid19_deaths_US.csv",'ReadVariableNames',1);
recoversUS=readtable("daily.csv",'ReadVariableNames',1);


confirmedUS(:,[1:4,8:11])=[];
deathsUS(:,[1:4,8:11])=[];



Mainland = confirmedUS.FIPS > 1000 & confirmedUS.FIPS < 56046;
uscounties_confirmed=confirmedUS(Mainland,:);

Mainland = deathsUS.FIPS > 1000 & deathsUS.FIPS < 56046;
uscounties_deaths=deathsUS(Mainland,:);

pop_table=uscounties_deaths(:,4);

% add state fips
state_fips=table(floor(uscounties_confirmed.FIPS/1000));
uscounties_confirmed=[state_fips uscounties_confirmed];
uscounties_deaths=[state_fips uscounties_deaths];



%% confirmed and deaths 

confirmed_case_data_array = table2array(uscounties_confirmed(:,5:end))';

deaths_case_data_array = table2array(uscounties_deaths(:,6:end))';

pop_array=table2array(pop_table)';

pop_array_multi=repmat(pop_array,length(confirmed_case_data_array(:,1)),1);

% normalize

norm_confirmed=confirmed_case_data_array./pop_array_multi;

norm_deaths=deaths_case_data_array./pop_array_multi;



%% Recovered if using real data (not applicable to states not reported)


daily_use=[recoversUS.date recoversUS.recovered,recoversUS.fips];
stat_fips=daily_use(:,3);

rows=stat_fips >0 & stat_fips<57;
recover1=daily_use(rows,:);
recover1(isnan(recover1))=0;

datefile=unique(recover1(:,1));
fips_file=unique(recover1(:,3));


M=zeros(length(datefile),length(fips_file));
for i=1:length(datefile)
       for j=1:length(fips_file)
ok=recover1(:,1)==datefile(i) & recover1(:,3)==fips_file(j);
A=recover1(ok,2);
if isempty(A)
     A=0;
end

M(i,j)=A;


        end
end

recover2= [fips_file';M];
zero_datefile=[0;datefile];
recover3=[zero_datefile recover2];
recover5=recover3';

new_structure=uscounties_confirmed(:,1:4);



for k=1:height(new_structure(:,1))
location=find(new_structure.Var1(k)==recover5(:,1));
 x=recover5(location,2:end);
 x_all(k,:)=x;
end

x_all_table=array2table(x_all);
recover_state=[new_structure x_all_table];


for n=1:length(recover5(:,1))
aaa=uscounties_confirmed.Var1==recover5(n,1);
bbb=uscounties_confirmed(aaa,5:end);
arraybbb=table2array(bbb);
sum_bbb=sum(arraybbb,1);
% bbb_all=repmat(sum_bbb,length(arraybbb(:,1)),1);
sum_bbb_all(n,:)=sum_bbb;

end


for k=1:height(new_structure(:,1))
location2=find(new_structure.Var1(k)==recover5(:,1));
 y=sum_bbb_all(location2,:);
 y_all(k,:)=y;
end

 y_all_table=array2table(y_all);
 confirmed_state=[new_structure y_all_table];
 
 confirmed_county_array=confirmed_case_data_array';
 confirmed_state_array=y_all(:,:);
recovered_state_array=x_all(:,:);

 
t=12;

[asdf,ghjk]=size(confirmed_state_array);
recovered_county_array=zeros(asdf,ghjk);

for p=t+1:ghjk
recovered_county_array(:,p)=recovered_state_array(:,p).*confirmed_county_array(:,p-t)./ confirmed_state_array(:,p-t);
end
 
recovered_county_array(isnan(recovered_county_array))=0;
recovered_case_data_array=recovered_county_array';
norm_recovered=recovered_county_array'./pop_array_multi;

%% Classify urban and rural



land_fips=readtable('land_by_FIPS.csv');
land_fips(1,:)=[];



fips=table2array(land_fips(:,2));
fips(length(fips)+1)=2158;
fips(length(fips)+1)=46102;
area=table2array(land_fips(:,3));
area(length(area)+1)=17081.43;
area(length(area)+1)=2093.9;
area(92)=446.4;
area(96)=2593.8;
area(97)=7762.1;
area(75)=12637;
area(258)=33.57;

fips_table=array2table(fips);
area_table=array2table(area);


known_fips=table2array(confirmed_state(:,2));




for i=1:length(known_fips)
[Lcob,Lia]=find(fips==known_fips(i));
area_all(i)=area(Lcob);
lcob_all(i)=Lcob;
end

known_pop=table2array(uscounties_deaths(:,5));

pop_density=known_pop./area_all';

rank=zeros(length(pop_density),1);


threshold=500;

for i=1:length(known_fips)
    if pop_density(i)>threshold
    rank(i)=1;
    elseif  pop_density(i)<=threshold && pop_density(i)>=0
      rank(i)=2;
    else
      rank(i)=3;
    end
end 

den_cat=[known_fips pop_density rank];


%% Save the data

csvwrite('../Data/COVID_outputs/normalized_confirmed_time_by_nodes.csv', norm_confirmed);
csvwrite('../Data/COVID_outputs/normalized_deaths_time_by_nodes.csv', norm_deaths);
csvwrite('../Data/COVID_outputs/normalized_recovery_time_by_nodes.csv', norm_recovered);
csvwrite('../Data/COVID_outputs/population_density_and_urbanClassification.csv', den_cat);






