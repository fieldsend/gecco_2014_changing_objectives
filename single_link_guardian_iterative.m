function [Y,Y_dom,I_a] = single_link_guardian_iterative(Y,Y_dom,I_a,n,chg,type_vector)

% [Y,Y_dom,I_a] = GECCO_single_link_guardian_iterative(Y,Y_dom,n,chg,sim)
%
% Code relates to:
% Fieldsend JE & Everson RM.
% Efficiently Identifying Pareto Solutions when Objective Values Change, 
% Genetic and Evolutionary Computation Conference, 
% GECCO'14, pages 605-612, Vancouver, Canada, 12th - 16th Jul 2014.
%
% Please cite the above work if you use this code
%
% Iteratively maintains a matrix of evaluated solutions (Y) as members are
% either changed (e.g. due to Maximum Likelihood location updates, dynamic 
% function change, etc) or new due to new locations being added. This is 
% done via simulating single link maintenance 
%
% Uses the dist2 function from the Netlab toolbox (freely available online)
%
% Inputs:
% 
% Y = matrix of maximum liklihood objective vectors (each row is a vector)
% Y_dom  = vector of single domination links
% I_a = indices of current elite (non-dominated) members of Y
% n = number of valid  members of Y
% chg = index of changed solution (if no solution has been changed, this
% should be 0, and the nth row of Y will be treated as a new solution)
% type_vector = vector of 6 elements holding the flag type for each of the
%      link maintenance checks. Please refer to the above paper for a 
%      definition of each option value. (Detailed in Algorithm 1 and 2.) 
%      Elements at indices 2 and 4-6 can take the values 1,2,3 or 4. 
%      Elements at indice 1 and 3 can take the values 1,2 or 3.
%
% Outputs:
%
% Y = matrix of maximum liklihood objective vectors (each row is a vector)
% Y_dom  = vector of single domination links 
% dom_comps = number of domination comparisons used in this iterative
%
% copyright Jonathan Fieldsend, 2013, 2014, 2016

if (chg>n)
   error('index of changed location cannot be higher than number of valid rows in Y'); 
end

[~,d] = size(Y);
% single iteration


% update state
if chg==0 % if no change has been made, then nth location is new
    [Y, Y_dom,I_a] = update_state_new_sample(Y, Y_dom, I_a, n, type_vector,d);
else % otherwse update due to change in objective vector of 'chg'th solution
    [Y, Y_dom,I_a] = update_state_resample(Y, Y_dom, I_a, n, chg, type_vector,d);
end
    
    

%------
function [X, X_dom,I_a] = update_state_resample(X, X_dom, I_a, n, chg, type_vector,d)

% function to update links if a single element at index 'chg' has been
% moved

if isempty(I_a)
    error('empty archive list');
end

% get indices of all solutions which have the changed member as a guide. 
% Uses Matlab's vector operations -- however in an OO language I would 
% suggest storing pointers in both directions rather than using a 'find' 
% call
% As dynamically extending/contracting vectors stored for each element in
% Matlab is not particularly efficient I've used find here, however I 
% should release a Java version at some point maintaining reference/array 
% of references for each guarded/guarding solution. This will negate the 
% need for a 'find' call

I_cc = find(X_dom(1:n)==chg);  
if sum(I_cc==chg)>0
    error('member seems to be guiding itself!');
end
    
ledge=0;
% if resampling a leading edge member
if sum(I_a==chg)==1
    I_a(I_a==chg)=[]; % remove from leading edge
    ledge=1; % flag to track removed from archive 
else
    old_dom = X_dom(chg); % can use previous dominator of changed to act 
                          % as guardian for those previously guided by chg
end

I_aa=I_a; % track indices after any removal
I_c = [I_a I_cc']; % those who new need a guide, plus archive estimate
I_ci = [I_cc' chg]; % those who need a new guide, plus changed
I_cii = [I_cc' chg X_dom(chg)]; % those who need a new guide, plus changed and changed previous guide
I_d = [I_cc' X_dom(chg)]; % those who need a new guide, plus changed previous guide

% stripping out any duplicates
I_ci = union(I_ci, I_a,'stable');
I_cii = union(I_cii, I_a,'stable');
I_d = union(I_d, I_a,'stable');

I_a=[]; % clear elite archive 

if ledge == 1 % was in E^t, so guided solutions may enter E^t+1
    r = zeros(length(I_c),1);
    k = zeros(length(I_c),1);
    % first see if any of the previously guided plus archive at t dominate 
    % the changed position
    for j=1:d
        r = r + (X(I_c,j)<=X(chg,j));
        k = k + (X(I_c,j)<X(chg,j));
    end
    z= sum(r==d & k>0);
    if (z==0) % chg not dominated
        I_a = [I_a chg]; % add to archive
        r = zeros(length(I_aa),1);
        for j=1:d
            r =r + (X(I_aa,j)>=X(chg,j));
        end
        % add back into archive non-dominated elements of I_aa
        I_a = [I_a I_aa(r~=d)];
        X_dom(I_aa(r==d)) = chg; % set chg as new guardian dominator for dominated elements of I_aa
    else % chg now dominated
        I_a = I_aa;
        II = find(r==d  & k>0);
        if type_vector(3) == 1
            X_dom(chg) = I_c(II(1));
        elseif type_vector(3)==2 % closest of guided by elite
            II = find(r==d);
            dd= dist2(X(I_c(II),:), X(chg,:));
            [~,ind2] = min(dd);
            X_dom(chg) = I_c(II(ind2));
        elseif type_vector(3)==3 % least used of guided by elite
            II = find(r==d);
            %II has indiced of all dominating
            if (length(II)==1)
                X_dom(chg) = I_c(II(1));
            else
                K = zeros(length(II),1);
                for i=1:length(K)
                    K(i) =  sum(X_dom(1:(n-1))==I_c(II(i)));
                end
                [~,ind2] = min(K);
                X_dom(chg) = I_c(II(ind2));
            end
        end
    end
    % now check to reassign guarded elements
    for k=1:length(I_cc)
        r = zeros(length(I_ci),1);
        for j=1:d
            r =r + (X(I_ci,j)<=X(I_cc(k),j));
        end
        z = sum(r==d);
        zi=find(r==d);
        if z==1 && zi<=length(I_cc)  % only doms itself
            I_a = [I_a I_cc(k)];
            X_dom(I_cc(k))=0; % mark not guided by anyone
        else
            r(k)=0; % can't be dominator for istelf
            II = find(r==d);
               
            if type_vector(4)==1
                if sum(X(chg,:)<= X(I_cc(k),:))~=d
                    X_dom(I_cc(k)) = I_ci(II(1));
                end
            elseif type_vector(4)==2
                X_dom(I_cc(k)) = I_ci(II(1));
            elseif type_vector(4)==3
                dd =dist2(X(I_ci(II),:), X(I_cc(k),:));
                [~,ind2] = min(dd);
                X_dom(I_cc(k)) = I_ci(II(ind2));
            elseif type_vector(4)==4
                %II has index of all dominating
                if (length(II)==1)
                    X_dom(I_cc(k)) = I_ci(II(1));
                else
                    K = zeros(length(II),1);
                    for i=1:length(K)
                        K(i) =  sum(X_dom(1:(n-1))==I_ci(II(i)));
                    end
                    [~,ind2] = min(K);
                    X_dom(I_cc(k)) = I_ci(II(ind2));
                end
            end
        end
    end
else % sampled point not from elite archive, so guarded solutions cannot end E^t+1
    r = zeros(length(I_aa),1);
    for j=1:d
        r =r + (X(I_aa,j)<=X(chg,j));
    end
    z= sum(r==d);
    if (z==0) % not dominated
        I_a = [I_a chg];
        r = zeros(length(I_aa),1);
        for j=1:d
            r =r + (X(I_aa,j)>=X(chg,j));
        end
        % add back into archive non-dominated elements of I_aa
        I_a = [I_a I_aa(r~=d)];
        X_dom(I_aa(r==d)) = chg;
        X_dom(chg)=0; % mark as non-domed by no longer having guide
    else % now dominated
        I_a = I_aa;
        II = find(r==d);
        if type_vector(5) ==1 % select first dominating of E
            if sum(X(X_dom(chg),:)<=X(chg,:))~=d
                X_dom(chg) = I_a(II(1));
            end
        else
            r = zeros(length(I_d),1);
            for j=1:d
                r =r + (X(I_d,j)<=X(chg,j));
            end
            II = find(r==d);
            
            if type_vector(5) == 2 % select first dominating of Y_chg U y* U E
                X_dom(chg) = I_d(II(1));
            elseif type_vector(5)==3 % closest of guided 
                if (length(II)>1)
                    dd= dist2(X(I_d(II),:), X(chg,:));
                    [~,ind2] = min(dd);
                    if (I_d(II(ind2))==chg)
                        error('indices should not be equal');
                    end
                    X_dom(chg) = I_d(II(ind2));
                else
                    X_dom(chg) = I_d(II(1));
                end
            elseif type_vector(5)==4 % least used of guided
                %II has indiced of all dominating
                if (length(II)==1)
                    X_dom(chg) = I_d(II(1));
                else
                    K = zeros(length(II),1);
                    for i=1:length(K)
                        K(i) =  sum(X_dom(1:(chg-1))==I_d(II(i)));
                    end
                    [~,ind2] = min(K);
                    X_dom(chg) = I_d(II(ind2));
                end
            end
        end
    end
    
    % now check to reassign guarded elements
    for k=1:length(I_cc)
        r = zeros(length(I_cii),1);
        for j=1:d
            r =r + (X(I_cii,j)<=X(I_cc(k),j));
        end
        z = sum(r==d);
        zi=find(r==d);
        if z==1 && zi<=length(I_cc)  % only doms itself
            I_a = [I_a I_cc(k)];
            X_dom(I_cc(k)) = 0; % mark as no longer having guide
        else
            r(k)=0; % can't be dominator for itself
            II = find(r==d);
            if type_vector(6)==1
                if sum(X(chg,:)<= X(I_cc(k),:))~=d
                    X_dom(I_cc(k)) = old_dom;
                end
            elseif type_vector(6)==2
                X_dom(I_cc(k)) = I_cii(II(1));
            elseif type_vector(6)==3
                dd =dist2(X(I_cii(II),:), X(I_cc(k),:));
                [~,ind2] = min(dd);
                X_dom(I_cc(k)) = I_cii(II(ind2));
            elseif type_vector(6)==4
                %II has indiced of all dominating
                if (length(II)==1)
                    X_dom(I_cc(k)) = I_cii(II(1));
                else
                    K = zeros(length(II),1);
                    for i=1:length(K)
                        K(i) =  sum(X_dom(1:(n-1))==I_cii(II(i)));
                    end
                    [~,ind2] = min(K);         
                    X_dom(I_cc(k)) = I_cii(II(ind2));
                end
            end
        end
    end
end

%------
function [X, X_dom,I_a] = update_state_new_sample(X, X_dom, I_a,n,type_vector,d)

% function to update links if a brand new location at index n has entered
% (members at 1 to n-1 already existing)

% compare to archive
r = zeros(length(I_a),1);
for j=1:d
    r =r + (X(I_a,j)<=X(n,j));
end

if sum(r==d)==0
    % if not dominated by elite set
    r = zeros(length(I_a),1);
    
    for j=1:d
        r =r + (X(I_a,j)>=X(n,j));
    end
    II = find(r==d);
    X_dom(I_a(II)) =  n;
    I_a(II) = [];
    I_a = [I_a n];
else
    % is dominated, so choose first elite member
    II = find(r==d);
    if type_vector(1) ==1 % first elite
        X_dom(n) = I_a(II(1));
        dom_index = I_a(II(1));
    elseif type_vector(1) ==2 % chose closest dominating
        dd= dist2(X(I_a(II),:), X(n,:));
        [~,ind2] = min(dd);
        dom_index = I_a(II(ind2));
        X_dom(n) = dom_index;
    elseif type_vector(1) ==3 %chose  dominating with fewest guided solutions
        K = zeros(length(II),1);
        for i=1:length(K)
            K(i) =  sum(X_dom(1:(n-1))==I_a(II(i)));
        end
        [~,ind2] = min(K);
        dom_index = I_a(II(ind2));
        X_dom(n) = dom_index;
    end
    if type_vector(2)>1
        I_c = find(X_dom(1:(n-1))==dom_index);
        I_c = [I_c' dom_index]; % indices of points that may do, and the elite that def does
        r = zeros(length(I_c),1);
        for j=1:d
            r =r + (X(I_c,j)<=X(n,j));
        end
        II = find(r==d);
        if type_vector(2)==2 % first of guided by elite
            X_dom(n) = I_c(II(1));
        elseif type_vector(2)==3 % closest of guided by elite
            dd= dist2(X(I_c(II),:), X(n,:));
            [~,ind2] = min(dd);
            X_dom(n) = I_c(II(ind2));
        elseif type_vector(2)==4 % least used of guided by elite
            %II has indiced of all dominating
            if (length(II)==1)
                X_dom(n) = I_c(II(1));
            else
                K = zeros(length(II),1);
                for i=1:length(K)
                    K(i) =  sum(X_dom(1:(n-1))==I_c(II(i)));
                end
                [~,ind2] = min(K);
                X_dom(n) = I_c(II(ind2));
            end
        end
    end
end


