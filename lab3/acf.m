function [ r_barlett, r_tukey] = acf( x )
% Estimates the acf of a signal x
% Assumes x is on the form (N,S) = (signal_length, number_of_signals)
% Returns both tukey and barlett estimates

N = size(x,1);
S = size(x,2);
r_tukey = zeros(2*N-1,S);
r_barlett = zeros(2*N-1,S);

for k = 1:N
    r = sum(x(1:N+1-k,:) .* x(k:end,:)) ;
    r_tukey(k,:) = r / (N + 1 - k);
    r_barlett(k,:) = r / N;
end

r_tukey(N+1:end,:) = flip(r_tukey(2:N,:), 1);
r_barlett(N+1:end,:) = flip(r_barlett(2:N,:), 1);


%r_tukey = [flip(r_tukey(2:end)), r_tukey];
%r_barlett = [flip(r_barlett(2:end)), r_barlett];

end



