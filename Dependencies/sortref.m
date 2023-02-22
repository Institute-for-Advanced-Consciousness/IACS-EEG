% Find the indicies of B in A (i.e., if 4.2 is the 9th element of A and the
% 1st element of B, returns 9 in the first element of output. Returns 
% vector C, the indicies needed to sort A according to B. Output from this 
% function can then be used to sort another vector (with entirely different 
% values) according to B (in combination with the Matlab sort() function).
% 
% Joel Frohlich, UCLA, 2021
%
% Examples
% 
% tmp = sort(dates); % dates for behavioral assessment
% NDX = sortref(dates,tmp); % find how the order changed when we sorted
% scores = scores(NDX); % sort corresponding behavioral scores according to dates 
% dates = tmp; % finish sorting dates
% 
% clearvars
% X = randperm(10); % vector to be sorted
% Y = randperm(10); % reference vector
% NDX = sortref(X,Y);
% all(X(NDX) == Y) % should return true

function[C] = sortref(A,B)

assert(size(A,1) == size(B,1),'Reference vector and vector to be sorted contain different number of values')
assert(size(A,2) == size(B,2),'Reference vector and vector to be sorted contain different number of values')

% transpose if column vector
if size(A,2) == 1
    A = A';
    B = B';
end

C = nan(size(A)); % allocation

for irow = 1:size(A,1)
    if iscell(A)
        if ~isempty(setdiff(A,B))
            warning('One vector contains values not contained in the other')
        end

        for icol = 1:size(A,2)
            C(irow,icol) = find(strcmp(A,B(irow,icol)));
        end
    else
        if any(unique(A) ~= unique(B))
            warning('One vector contains values not contained in the other')
        end

        for icol = 1:size(A,2)
            if isnan(B(irow,icol))
                try 
                    C(irow,icol) = find(isnan(A(irow,:)));
                catch
                    if length(find(isnan(A(irow,:)))) > 1 % if more than one match
                        matches = find(isnan(A(irow,:)));
                        for j = 1:length(matches)
                            if ~ismember(matches(j),C(irow,:)) % if this index hasn't been added yet ...
                                C(irow,icol) = matches(j); % ... add it to list
                                break
                            end
                        end
                    end
                    if isnan(C(irow,icol))
                        error('Failed to assign index')
                    end
                end  
            else
                try
                    C(irow,icol) = find(A(irow,:) == B(irow,icol));
                catch
                    if length(find(A(irow,:) == B(irow,icol))) > 1 % if more than one match
                        matches = find(A(irow,:) == B(irow,icol));
                        for j = 1:length(matches)
                            if ~ismember(matches(j),C(irow,:)) % if this index hasn't been added yet ...
                                C(irow,icol) = matches(j); % ... add it to list
                                break
                            end
                        end
                    end
                    if isnan(C(irow,icol))
                        error('Failed to assign index')
                    end
                end
            end
        end
    end
end

end






