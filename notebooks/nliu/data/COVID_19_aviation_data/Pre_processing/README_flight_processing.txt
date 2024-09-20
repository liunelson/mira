Before running "Flight_preprocessing.m" , you need to have the following in "Transportation_air_network_dataset" folder :

1. "Jan200.csv", "Feb200.csv".... "Jul200.csv" (we will keep updating). 
2. "Fleet.csv" which is a dataset download from BTS with all registered aircraft with configurations(seat number)

3."Airport_fips_northeastern.csv" is the fips Code of the airport in the northeastern area.

4. Run COVID_processing first to ensure success.

You will output the following in "Data" folder:

1. "A3_Profile.mat" is the flight adjacency matrix in northeastern area between counties.
The number in the matrix is the number of seats between two city pairs/ the maximum number of seats between in all city pairs 

2. "List_of_flight.csv" is the day-by-day flight information (since Jan 22, 20 to Jul 31, 20) between any city pairs in northeastern area.



 

