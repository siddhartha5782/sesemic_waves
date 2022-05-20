clear; close all;
initstate(2999);

NN = 201 ;
L = 1000 ;

garch = zeros(L, NN);

for ll = 1:NN

    e = randn(L, 1);

	    %% garchsimulate is from UCSD garch simulation code
    [noise, H] = garchsimulate(length(e),[1, 0.2, 0.2, 0.3],1,2);
    garch(:, ll) = noise;

end

save garch_noise garch


