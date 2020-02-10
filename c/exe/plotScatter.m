%plot the rate and eq probability scatter data
%can also plot pareto data and 1 interaction type data on top
clear; clf;

%set flags for which data to include in plot
pareto = 0;
intCompare = 1;

%get the rate and probability data
data = importdata("rateEqScatter.txt"); 
prob = data(:,1); rate = data(:,2); 

%plot the basic data first
hold on
scatter(prob,rate,'b')
xlabel("Equilibrium Probability");
ylabel("Average Transition Rate");

%check if using pareto and if the file exists
if pareto == 1
    if exist("rateEqPareto.txt", "file")
        dataP = importdata("rateEqPareto.txt");
        p1 = dataP(:,1); r1 = dataP(:,2);
        scatter(p1,r1,'r')
    else
        error("The pareto file does not exist");
    end
end

%check if using intCOmpare and if file exists. file name needs to be
%manually input
if intCompare == 1
    if exist("triangleScatter.txt", "file")
        dataC = importdata("triangleScatter.txt");
        p2 = dataC(:,1); r2 = dataC(:,2);
        scatter(p2,r2,'g')
    else
        error("The comparison file does not exist");
    end
end

hold off




