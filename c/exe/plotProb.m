%plot the hitting probabilities
clear;

d = importdata("hittingProb.txt"); 

N = d.data;
T = d.textdata;

num_states = length(T)-1;

kappa = N(1,:);
M = length(kappa);
eq1 = [48.5,39.9,11.4]/100; color1 = ["blue", "green", "red"];

%for 6 particles
figure(1)
hold on
for i=1:num_states
    plot(kappa,N(i+1,:),"Color",color1(i))
end
for i=1:num_states
    plot(kappa,ones(M,1)*eq1(i),"--","Color",color1(i))
end
hold off
xlabel("Sticky Parameter");
ylabel("Hitting Probability");
legend("Trapezoid","Chevron","Triangle")

figure(2)
kappa_final = 10;
hold on
for i=1:num_states
    plot(kappa(1:M/2+kappa_final-10),N(i+1,1:M/2+kappa_final-10),"Color",color1(i))
end
for i=1:num_states
    plot(kappa(1:M/2+kappa_final-10),ones(M/2+kappa_final-10,1)*eq1(i),"--","Color",color1(i))
end
xlabel("Sticky Parameter");
ylabel("Hitting Probability");
hold off
legend("Trapezoid","Chevron","Triangle")
