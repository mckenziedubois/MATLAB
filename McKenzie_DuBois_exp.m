function McKenzie_DuBois_exp(starting_number,growth_rate,gen_plot)

sim_growth(starting_number,growth_rate,gen_plot,10) %called our function sim_grwoth and varied the value of n_points
sim_growth(starting_number,growth_rate,gen_plot,100)
sim_growth(starting_number,growth_rate,gen_plot,1000)
sim_growth(starting_number,growth_rate,gen_plot,10000)

nine_doubles(growth_rate) %called our function nine_doubles

end


function [sampling_times,cell_count]= sim_growth(starting_number,growth_rate,gen_plot,n_points)

estimate_six_doubles = 6*(69.31471806/growth_rate); %found the amount of time it takes to complete six doublings
times = linspace(0,estimate_six_doubles,n_points); %created an array of evenly spaced values from 0 to the time above

%npoints defines how many values we have between 0 and estimate_six_doubles

count = zeros(length(times)); %an array of zeros, with a zero for every value of time 
count(1) = starting_number; 

    for i = 2:length(times) %initiate a for loop that loops through all the values of time
        dt = times(i) - times(i-1); %using each value of time and the previous value of time to find the dt
        count(i) = count(i-1) + (count(i-1) * dt * (growth_rate/100)); %uses this delta t to find the count at the particular value of t 
    end
    
    if gen_plot == 1 % generates the plot if asked for 
        figure(1);
        hold on
        plot(times,count);
        title('Cell Count');
        xlabel('Time (hours)')
        ylabel('Cell Count')
        legend({sprintf('Growth rate = %0.0f%', growth_rate), sprintf('Starting number of cells = %f%', count)})
    end
  
sampling_times = times(length(times)); %grabs the final value of time 
cell_count = count(length(count)); %grabs the final count value

[sampling_times,cell_count] % returns the two values listed in line 38 and 39
end
    

function nine_doubles(growth_rate)

estimate_nine_doubles = 9*(69.31471806/growth_rate);  %found the amount of time it takes to complete nine doublings

starting_values = zeros(10000,1); %created an array of evenly spaced values from 0 to 10000 
ending_values = zeros(length(starting_values),1); % created an array of zeros for the ending values 
    
for x = 1:length(starting_values) % initiated a for loop to run 10000 times assigning each zero in starting values a value of 1 to 10000
    starting_values(x) = x;
end

for indx = 1:length(starting_values) % initiated a for loop to run 10000 times  
    ending_values(indx) = starting_values(indx)*exp((2/100)*estimate_nine_doubles); % calculated the ending value at each starting value 
end

figure(2); 
plot(starting_values,ending_values); % plotted a figure of starting values versus ending values 
title(sprintf('# of cells after 9 doublings v. initial cells'));
    xlabel('initial # cells')
    ylabel('# of cells after 9 doublings')
    b = polyfit(starting_values,ending_values, 1); % found the slope of the line 
    slope = b(1);
    legend("Slope:" + slope); %aded the slope to the legend 

end 

