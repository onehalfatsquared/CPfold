%%run the mfpt estimator - in parallel. recieves the numerator and
%%denominator of the estimator as well as the non-normalized row of the embeded
%%transition matrix. stores these in a cell 
% storage - (1,2,3,4) = (num, denom, transmatx, times ran)
clear

%set parameters
N = 7;
num_test = 1;

%import the graph structure with all states
%Check if database exists for this value of N
filename = strcat('N', num2str(N), 'DB.mat');
if exist(filename, 'file')
    %file exists, load it
    load(filename);
else
    error('Requested database does not exist. Initialize it first');
end

%import the mfpt cell storage
%Check if database exists for this value of N
filename = strcat('N', num2str(N), 'mfpt.mat');
if exist(filename, 'file')
    %file exists, load it
    load(filename);
else
    %initialize the database
    L = length(DB);
    Pempty = zeros(1,L);
    for i = 1:L
        mfpt{i,1} = 0; mfpt{i,2} = 0; mfpt{i,4} = 0; mfpt{i,3} = Pempty;
    end
end

L = length(mfpt);
err
%pick a state that hasnt been run yet. do test
for i = 1:L
    if mfpt{i,4} < num_test
        num = mfpt{i,1}; den = mfpt{i,2}; P = mfpt{i,3};
        [num, den, P] = estimateMFPTsingle(N, i, num, den, P);
        mfpt{i,1} = num; mfpt{i,2} = den; mfpt{i,3} = P; mfpt{i,4} = mfpt{i,4} + 1;
    end
    i %output progress
end

%save('N7mfpt.mat', 'mfpt')