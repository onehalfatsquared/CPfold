%tries to extract more data from the scatter plot files
%ex: get the values from varying AB at fixed AA and BB
clear; clf;

%get the rate and probability data
data = importdata("rateEqScatter.txt"); 

%set the kappa values used to generate the scatter data (may change)
%comes from a logarithmic scale
base = 0.5;
multiplier = 1.6;
M = 33;
indices = 0:M-1;
kappa = 0.5 * multiplier.^indices;

%set the test case
test = 2;


%set a sequence to extract 
%test 1 - fixing AA and AB
if (test == 1)
    AAfix = 1;
    ABfix = 15;
    BB = (1:M) + ((AAfix-1) * M^2 + (ABfix-1) * M);
    probs = data(BB,1); rates = data(BB,2);
end

%test 2 - fixing AA and BB
if (test == 2)
    AAfix = 1;
    BBfix = 33;
    AB = ((0:M-1) * M) + ((AAfix-1) * M^2 + BBfix);
    probs = data(AB,1); rates = data(AB,2);
end

%test 3 - fixing AB and BB
if (test == 3)
    ABfix = 10;
    BBfix = 10;
    AA = ((0:M-1) * M^2) + ((ABfix-1) * M + BBfix);
    probs = data(AA,1); rates = data(AA,2);
end

%plot results
figure(1)
plot(log(kappa), probs)
xlabel("log(\kappa)")
ylabel("Equilibrium Probability")
figure(2)
plot(log(kappa), rates)
xlabel("log(\kappa)");
ylabel("Average Transition Rate")






