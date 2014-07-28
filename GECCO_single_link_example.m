function [Observed,Truth,Y_dom,I_a,I_n,t]=GECCO_single_link_example(n,d,type_vector,S,seed)

% function shows example single link generation and tracking
%
% alternate time steps a location is generated or moved
%
% Code relates to:
% Fieldsend JE & Everson RM.
% Efficiently Identifying Pareto Solutions when Objective Values Change, 
% Genetic and Evolutionary Computation Conference, 
% GECCO'14, pages 605-612, Vancouver, Canada, 12th - 16th Jul 2014.
%
% Please cite the above work if you use this code
%
% Inputs
%
% n = number of time steps - must be even, otherwise rounded to lower even
%      number
% d = number of objectives
% type_vector = vector of 6 elements holding the flag type for each of the
%      link maintenance checks. Please refer to the above paper for a 
%      definition of each option value. (Detailed in Algorithm 1 and 2.) 
%      Elements at indices 2 and 4-6 can take the values 1,2,3 or 4. 
%      Elements at indice 1 and 3 can take the values 1,2 or 3.
% S = simulation number (valid values are 1, 2 3 and 4). Again, definitions
%      are provided in the paper
% seed = random number generator seed
%
% Outputs:
%
% Observed = matrix of maximum liklihood objective vectors (each row is a vector)
% Truth = matrix of actual (noise-free) objective vectors
% Y_dom  = vector of single domination links 
% dom_comps = number of domination comparisons used in this iterative
% t = time steps taken
%
% copyright Jonathan Fieldsend, 2013,2014


rng(seed);
n2 = floor(n/2);
t=1;

display('starting');

% preallocate matrices for speed
Observed = zeros(n2,d); % estimated objective values
Truth = zeros(n2,d); % actual objective values
Y_dom = zeros(n2,1); % index of dominating element which 'guards' each row vector (0 if not dominated)
I_n = ones(n2,1); % number of resamples per location

display('matrices preallocated');

% make initial
Truth(1,:) = zeros(1,d);
Observed(1,:) = Truth(1,:) + randn(1,d)*0.1; % observed is truth plus noise term
 
I_a = t; % with a single instance evaluated, this is also the Pareto estimate

for i=1:n2-1
    % change one
    [Observed,Truth, Y_dom, I_a, I_n, t] = change_a_member(Observed, ...
        Truth, Y_dom, I_a, I_n, S, type_vector, t, d, i);
    if rem(t,1000)==0
        fprintf('Time steps: %d\n',t);
    end
    % add one
    [Observed,Truth, Y_dom, I_a, I_n, t] = add_a_member(Observed, ...
        Truth, Y_dom, I_a, I_n, S, type_vector, t, d, i);
end

[Observed,Truth, Y_dom, I_a, I_n, t] = change_a_member(Observed, ...
    Truth, Y_dom, I_a, I_n, S, type_vector, t, d, n2);
fprintf('Final time step: %d\n',t);

%----------

function [Observed,Truth, Y_dom, I_a, I_n, t] = change_a_member(Observed, ...
    Truth, Y_dom, I_a, I_n, S, type_vector, t, d, i)

if S <3
    chg = randperm(i);
    chg = chg(1); % randomly choose a member to change
else
    chg = randperm(length(I_a));
    chg = I_a(chg(1)); % randomly choose an elite member to change
end
I_n(chg) = I_n(chg) + 1;
% observed is updated with another sample of the truth plus noise term
Observed(chg,:) = Observed(chg,:) + ((Truth(chg,:) + randn(1,d)*.02 - Observed(chg,:))/I_n(chg)); % incremental mean estimate
[Observed,Y_dom,I_a] = single_link_guardian_iterative(Observed,Y_dom,I_a,i,chg,type_vector);
t= t+1;



%----------

function [Observed,Truth, Y_dom, I_a, I_n, t] = add_a_member(Observed, ...
    Truth, Y_dom, I_a, I_n, S, type_vector, t, d, i)

if S==1 || S==3
    Truth(i+1,:) = randn(1,d); %new solution randomly generated
else
    opt = randperm(length(I_a));
    opt = I_a(opt(1)); % select a current optimal estimate at random
    Truth(i+1,:) = Truth(opt,:) + randn(1,d)*0.25; % new solution generated in vicinity of current estimate
end
Observed(i+1,:) = Truth(i+1,:) + randn(1,d)*0.1; % observed is truth plus noise term

[Observed,Y_dom,I_a] = single_link_guardian_iterative(Observed,Y_dom,I_a,i+1,0,type_vector);
t= t+1;