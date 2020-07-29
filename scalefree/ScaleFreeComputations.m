%Computations associated with scale free networks
%Author: Theju Jacob

clc;
clear all;
% close all;

lambda = 4.0; minconn = 7; maxconn = 117;
C = (lambda - 1)*minconn^(lambda - 1);
scalefreeprob = zeros(1,maxconn - minconn + 1);%[ 0.9 0.346989 0.145355 0.0652843 0.0311145 0.0156074 0.00818545 0.0044643 0.00252068 0.00146793 0.000878906 0.000539572 0.000338856 0.000217252 0.000141948 9.43718e-05 6.37542e-05 4.37124e-05 3.03852e-05 2.13925e-05 1.52416e-05 1.09806e-05 7.99361e-06 5.87629e-06 4.35967e-06 3.26259e-06 2.4616e-06 1.87165e-06 1.43353e-06 1.1056e-06 8.58307e-07 6.70508e-07 5.26926e-07 4.16445e-07 3.30914e-07 2.64312e-07 2.1216e-07 1.71105e-07 1.38621e-07 1.12793e-07 9.216e-08 7.56033e-08 6.226e-08 5.14617e-08 4.26879e-08 3.55317e-08 2.96731e-08 2.48596e-08 2.08912e-08 1.76085e-08 1.48844e-08 1.26167e-08 1.07232e-08 9.1377e-09 7.80626e-09 6.68512e-09 5.73856e-09 4.93736e-09 4.25749e-09 3.67918e-09 3.18612e-09 2.76477e-09 2.40391e-09 2.09418e-09 1.82779e-09 1.5982e-09 1.39993e-09 1.22839e-09 1.07968e-09 9.50543e-10 8.3819e-10 7.40274e-10 6.54793e-10 5.80045e-10 5.14576e-10 4.57144e-10 4.06685e-10 3.62285e-10 3.23159e-10 2.88631e-10 2.58117e-10 2.31115e-10 2.07188e-10 1.85957e-10 1.67095e-10 1.50316e-10 1.35372e-10 1.22046e-10 1.10149e-10 9.95155e-11 9e-11 8.14758e-11 7.38313e-11 6.69685e-11 6.08008e-11 5.52522e-11 5.02555e-11 4.57514e-11 4.16874e-11 3.8017e-11 3.46989e-11];;
j = minconn;
for i = 1:(maxconn-minconn+1)
    scalefreeprob(i) = C*1/(j^lambda)
    j = j+1;
end
%scalefreeprob = [0.3 0.190181 0.12936 0.0926449 0.0690144 0.053033 0.0417892 0.0336196 0.0275225 0.0228679 0.019245 0.0163775 0.0140742 0.0122001 0.0106577 0.009375 0.00829847 0.00738737 0.00661038 0.00594317 0.00536656 0.00486534 0.00442728 0.00404251 0.00370298 0.00340207 0.00313431 0.00289515 0.00268078 0.00248799 0.00231407 0.0021567 0.00201392 0.00188403 0.00176557 0.00165728 0.00155807 0.00146698 0.00138317 0.00130591 0.00123457 0.00116856 0.00110739 0.00105061 0.000997829 0.000948683 0.000902861 0.000860078 0.000820081 0.00078264 0.000747549 0.000714622 0.00068369 0.0006546 0.000627215 0.000601407 0.000577061 0.000554073 0.000532347 0.000511795 0.000492337 0.0004739 0.000456414 0.000439819 0.000424056 0.000409073 0.000394821 0.000381255 0.000368332 0.000356014 0.000344265 0.000333052 0.000322344 0.000312111 0.000302328 0.000292969 0.00028401 0.00027543 0.000267209 0.000259327 0.000251767 0.000244512 0.000237546 0.000230855 0.000224425 0.000218243 0.000212296 0.000206574 0.000201066 0.000195761 0.00019065 0.000185724 0.000180974 0.000176393 0.000171972 0.000167705 0.000163585 0.000159605 0.000155759 0.000152042 0.000148448];
figure;axes('FontSize',24);
plot(minconn:maxconn,scalefreeprob);grid on; title('Connections vs Probability, Scale Free');xlabel('Number of connections');ylabel('Probability of connections');

% actualconn = [28 14 42 14 119 7 70 42 49 21 14 7 28 245 14 14 28 7 21 7 42 133 7 91 7 7 77 35 91 56 35 7 7 7 7 21 14 7 28 14 14 14 14 63 7 28 21 42 14 28 14 7 28 35 49 84 14 7 98 7 7 28 14 14 35 14 14 77 49 7 7 105 7 14 14 14 77 98 7 28 14 56 7 7 21 182 7 14 7 7 7 49 7 7 49 21 42 21 14 77 7 14 28 21 14 21 35 7 7 21 14 7 7 7 14 7 7 14 56 98 56 112 28 7 14 42 7 14 7 525 14 7 28 7 21 7 98 21 63 7 7 14 350 147 49 21 7 14 63 56 14 7 63 7 14 14 7 21 21 14 28 7 7 14 7 63 84 7 77 14 133 7 28 91 21 21 28 49 28 14 21 98 35 7 42 14 14 42 105 21 357 42 14 21 21 14 42 35 182 7 7 14 7 14 77 14 21 28 35 105 42 21 28 49 28 7 14 28 7 63 28 21 28 7 7 28 28 28 14 21 35 154 14 7 42 7 7 56 119 49 7 14 7 14 28 49 7 7 35 42 14 7 14 7 14 168 7 14 21 56 14 7 14 42 56 112 14 21 7 63 14 49 7 7 7 14 196 175 49 21 105 21 14 7 28 21 7 28 7 49 7 35 14 7 21 14 14 21 35 28 28 14 21 77 49 28 7 7 7 14 28 14 14 35 56 14 14 21 245 21 14 21 21 42 35 35 35 28 28 112 21 7 7 126 42 7 7 133 14 112 14 35 21 14 42 112 217 28 14 119 21 21 49 56 133 7 14 28 14 21 7 35 21 42 14 14 7 21 28 462 7 7 77 140 7 63 14 14 21 28 7 35 7 63 35 14 21 7 7 21 77 189 56 7 7 7 7 35 112 14 294 70 21 14 14 7 7 21 7 7 7 49 693 21 7 14 7 7 56 7 21 21 77 35 28 14 7 70 7 21 28 70 35 21 14 105 7 7 98 42 7 56 28 28 14 7 14 56 14 7 133 28 7 63 14 42 7 574 63 21 112 7 14 35 28 7 112 35 42 7 14 14 7 28 21 14 7 14 21 35 14 21 28 7 14 7 28 7 21 14 7 7 7 21 14 21 42 56 14 28 28 203 42 28 28 14 42 21 28 7 147 133 14 7 14 14 56 21 203 7 14 14 7 35 7 63 14 70 7 14 14 14 28 7 14 70 21 21 259 7 7 28 28 14 14 7 7 42 14 14 7 70 56 14 42 14 7 21 7 56 7 7 182 21 63 7 28 7 63 7 14 14 49 7 98 21 70 7 42 21 21 49 7 28 70 21 7 7 126 35 7 21 49 7 14 21 21 28 14 70 7 7 14 14 7 7 14 21 91 63 105 203 7 35 7 28 7 14 7 7 7 98 14 21 35 14 28 7 42 63 7 14 7 7 14 70 7 7 7 21 28 7 14 7 7 7 14 42 112 14 21 21 7 21 28 7 42 28 98 28 56 28 42 91 63 84 14 7 49 35 14 28 56 21 49 42 7 63 63 21 7 7 28 42 21 21 42 7 63 14 77 28 14 28 77 49 7 28 28 28 28 7 21 28 63 28 14 35 28 42 63 14 35 28 7 14 35 63 21 14 70 7 7 7 7 14 7 7 14 7 35 14 35 14 7 14 28 35 105 7 203 14 7 77 42 7 21 14 168 7 28 14 21 56 35 7 28 14 77 14 21 7 322 84 133 21 7 28 7 63 70 7 28 7 140 56 28 91 28 42 56 70 7 35 98 21 14 14 168 42 119 7 21 7 7 14 7 35 21 196 21 14 126 21 126 7 84 7 105 63 14 7 42 7 14 7 21 14 42 7 14 98 56 7 28 7 7 28 42 21 7 21 21 7 168 14 14 28 14 35 28 21 28 21 28 7 70 35 7 14 7 98 56 42 7 7 28 84 28 7 21 42 98 14 7 7 14 35 7 42 7 7 21 14 14 21 7 14 21 21 56 35 91 7 21 14 42 28 14 28 7 21 35 91 14 49 14 21 7 21 49 7 119 21 14 7 7 7 7 42 14 7 56 7 14 210 7 7 7 14 49 21 14 42 7 14 7 7 28 21 147 7 42 42 7 56 7 70 112 21 7 7 14 7 28 7 49 14 7 98 21 14 35 7 7 98 35 63 21 56 21 35 119 21 7 14 7 14 14 21 14 21 7 42 28 7 14 7 42 21 28 7 49 28 14 56 98 28 77 14 70 42 7 21 7 35 56 77 7 28 21 7 168 63 42 42 42 112 28 7 7 70 14 7 21 56 49 14 14 63 7 7 21 217 7 7 7 35 14 105 21 49 21 14 42 168 7 315 7 14 14 49 28 21 14 70 7 14 14 7 35 21 14 42 7 63 21 49 42 21 14 21 7 7 28 14 7 35 14 7 56 7 7 7 56 14 21 63 7 175 7 7 140 35 7 161 259 28 42 7 14 14 7 7 7 7 21 7 77 7 14 14 7 63 35 14 63 14 28 7 7 14 14 7 56 14 7 21 7 42 63 14 91 21 21 14 21 7 14 42 14 70 35 14 7 21 28 42 7 14 105 91 7 21 7 14 14 49 7 189 21 21 7 42 7 7 105 28 35 28 56 7 42 7 14 77 14 70 84 7 14 21 7 14 7 77 7 14 21 7 28 14 7 154 49 14 42 21 35 217 56 56 105 182 35 35 63 105 28 7 7 21 56 35 28 7 28 343 28 28 14 49 21 42 7 28 35 14 42 7 42 7 7 49 7 49 14 42 35 35 7 35 14 63 28 35 7 28 14 21 70 21 133 7 49 42 21 7 21 42 42 70 14 14 7 14 42 7 21 42 21 7 14 14 14 35 21 49 14 7 28 56 35 7 378 14 7 28 14 14 14 189 63 42 7 14 21 7 21 56 28 7 14 49 28 70 21 7 7 28 7 21 7 7 7 28 14 98 35 119 7 28 35 28 56 14 28 42 7 7 14 14 7 21 7 21 7 7 35 7 98 7 21 7 14 28 7 7 49 21 42 42 28 35 7 14 7 112 14 70 7 56 14 42 168 7 14 98 49 35 21 21 28 21 21 35 63 7 14 21 14 98 42 21 42 7 91 35 28 14 21 14 28 98 14 28 175 14 7 77 21 21 14 7 7 14 63 7 49 7 28 7 14 21 7 7 21 49 21 7 7 7 21 21 35 28 14 63 49 70 7 21 28 14 70 28 35 63 21 14 35 7 7 21 98 14 7 7 21 7 21 7 28 28 35 21 140 14 7 42 70 21 7 28 63 42 28 112 7 7 14 35 7 28 14 7 49 42 14 14 91 21 35 21 21 84 35 14 14 210 28 7 7 7 7 49 7 7 147 7 7 7 7 70 35 14 14 14 21 21 77 14 7 21 77 14 35 7 84 21 28 14 14 7 7 14 56 70 84 7 21 91 63 21 42 28 56 21 14 70 14 21 7 231 21 35 42 28 7 14 21 14 56 63 63 56 7 7 42 147 14 14 7 21 7 7 7 28 28 7 14 455 28 21 21 7 7 49 14 14 7 28 7 14 133 14 14 7 7 70 14 777 7 56 7 21 7 28 7 7 112 63 308 7 7 7 56 91 14 49 91 35 14 21 7 35 21 28 7 7 21 301 56 14 7 112 28 21 14 28 14 70 35 14 49 21 14 42 182 49 28 14 14 7 70 21 7 84 49 14 35 21 7 21 14 7 28 28 14 7 14 14 7 7 21 49 7 63 21 7 7 14 35 28 28 14 7 7 7 21 7 7 21 7 7 7 35 14 14 14 14 21 21 28 28 224 7 28 21 35 77 14 14 7 7 35 7 70 21 28 154 63 14 21 28 14 7 7 21 21 21 49 7 7 7 14 7 7 7 42 7 28 14 14 21 14 28 21 63 7 287 7 14 7 14 119 28 7 49 14 21 21 7 105 7 7 14 7 14 14 98 35 7 7 21 7 49 21 28 7 42 70 35 77 21 7 7 21 35 21 14 35 28 112 14 14 21 7 28 7 21 14 42 7 35 35 77 126 7 21 35 7 49 21 147 7 21 21 63 42 7 28 7 7 14 28 63 14 98 84 7 14 42 7 126 77 21 77 119 42 14 28 7 7 14 56 7 7 35 21 105 28 14 42 147 21 7 196 14 42 77 56 35 105 14 21 42 7 7 14 91 14 14 7 21 28 7 56 7 7 35 49 28 7 14 42 7 56 21 7 28 7 21 77 14 21 63 7 14 14 35 14 7 112 7 49 84 42 7 77 112 84 28 21 7 210 7 14 42 63 7 14 21 14 7 14 77 35 7 7 21 14 14 168 42 7 21 35 77 7 42 77 35 63 308 7 49 14 7 28 49 7 28 7 14 14 14 28 14 56 28 112 7 126 7 14 14 14 70 42 7 14 35 21 112 77 84 21 56 35 21 28 42 35 7 7 49 7 14 7 28 98 28 7 14 217 7 7 119 21 217 21 7 21 14 7 21 56 70 84 392 7 7 28 42 21 7 35 7 14 56 7 42 14 14 28 49 28 14 7 63 14 35 21 7 14 28 28 42 7 21 7 175 14 7 14 35 7 14 161 21 28 91 14 21 14 28 14 63 434 21 28 21 35 7 525 637 7 7 7 28 14 28 98 14 7 7 7 14 14 7 35 77 7 28 21 35 7 14 14 35 14 35 7 7 35 7 49 7 7 56 7 21 21 14 777 21 70 98 77 7 63 35 49 49 7 7 7 28 49 28 56 35 7 14 14 21 14 14 49 7 42 14 98 35 14 182 70 84 49 21 42 28 21 14 7 7 224 91 35 7 35 35 21 21 28 63 98 7 7 14 7 21 21 7 28 350 7 28 56 14 7 14 7 7 14 28 7 98 175 49 77 7 224 7 7 63 7 28 126 7 77 7 7 182 49 49 21 35 14 7 7 56 42 28 7 14 7 7 35 14 35 21 14 70 28 7 14 49 112 7 28 14 14 84 42 63 21 42 35 42 21 84 49 105 14 14 14 28 14 42 7 7 189 14 14 35 21 28 126 98 28 7 70 7 49 42 35 7 7 21 21 7 49 56 56 14 14 7 42 7 112 28 35 35 182 133 14 28 21 259 49 21 49 35 14 112 7 42 7 7 7 56 21 14 112 63 14 14 14 77 98 21 42 14 28 7 14 42 7 7 70 21 14 777 21 7 63 56 14 14 7 154 35 77 14 14 14 63 56 35 77 7 14 14 21 161 21 28 7 42 21 7 42 21 7 14 42 14 28 7 14 91 49 35 7 49 98 7 7 7 14 7 7 21 126 28 7 21 98 14 7 42 28 21 49 21 35 21 70 98 49 7 7 28 21 21 28 28 7 147 21 28 14 7 7 98 14 14 14 7 7 70 28 7 28 14 14 14 21 7 7 14 28 7 14 147 7 49 7 21 42 56 21 196 77 84 35 7 35 28 21 42 14 7 28 7 56 14 14 14 14 14 7 49 98 777 28 28 49 7 56 7 21 7 7 21 14 7 70 21 14 28 35 98 35 14 119 63 147 7 28 7 35 14 14 7 28 14 42 56 42 21 7 7 14 105 7 14 7 28 14 21 14 56 7 49 84 21 14 21 14 7 84 21 70 14 21 21 14 14 14 14 14 14 7 63 28 7 77 7 14 14 7 21 7 14 21 21 777 21 42 14 42 14 49 35 35 14 21 14 21 7 49 7 14 21 42 56 63 14 70 14 7 28 7 28 7 7 49 63 42 105 14 70 42 14 21 49 147 49 14 42 7 35 126 84 7 7 14 77 91 14 14 21 7 7 14 28 42 105 28 14 7 28 77 35 42 42 77 70 7 7 7 42 35 14 133 7 14 119 14 21 70 7 7 70 21 35 49 14 7 35 7 21 14 14 7 413 14 21 42 7 35 28 28 21 7 21 7 21 7 35 21 133 7 7 14 14 21 7 56 14 28 70 56 21 98 14 35 14 42 28 14 7 42 7 7 14 14 21 7 21 14 7 35 7 105 14 7 35 7 7 42 7 56 7 21 28 98 35 7 196 7 14 7 7 42 42 14 21 56 42 287 7 7 7 14 28 7 14 7 161 7 63 7 7 28 7 7 42 7 7 35 91 7 7 49 7 28 56 49 14 56 7 7 56 14 28 21 98 35 14 7 35 63 7 28 7 21 14 14 21 7 49 7 7 21 21 14 7 7 28 14 28 42 28 42 21 28 35 77 7 56 14 21 42 14 7 7 133 14 14 14 7 70 35 14 7 28 42 7 7 28 7 14 21 14 70 7 21 14 7 14 7 14 77 21 7 7 14 14 21 14 14 7 14 77 35 70 14 63 7 7 77 14 133 7 7 7 7 14 7 21 14 42 35 14 14 7 7 7 21 14 7 14 63 14 14 21 28 21 35 7 7 70 7 35 21 35 14 7 42 14 49 28 14 56 21 7 7 7 14 21 217 14 7 14 21 7 7 140 224 14 21 14 21 14 70 7 21 14 21 28 49 14 14 7 77 21 112 7 28 70 63 7 7 203 42 28 7 7 7 42 7 21 91 14 98 7 7 42 175 7 21 7 42 7 35 49 7 119 91 7 28 42 14 14 21 42 49 21 77 7 7 21 7 35 119 21 7 14 7 7 28 28 7 91 7 21 28 28 28 14 28 42 28 7 14 21 28 28 7 42 7 49 7 14 14 70 35 28 63 28 7 14 42 7 105 28 21 28 35 49 14 35 112 21 7 7 35 7 21 7 21 14 49 7 70 7 14 21 21 7 238 7 56 28 7 84 7 77 35 91 7 28 35 28 14 63 42 7 91 63 14 42 175 7 7 14 7 7 14 7 21 21 56 7 77 7 28 7 14 49 42 147 7 14 14 7 49 21 14 7 7 21 49 21 7 14 7 35 28 14 14 70 35 49 63 7 70 21 70 21 28 21 28 7 28 7 63 14 7 7 28 7 7 14 7 42 14 14 49 77 21 14 7 28 56 112 42 7 49 7 63 49 7 301 35 35 49 42 7 7 7 42 7 21 175 14 7 21 14 28 7 21 21 7 42 14 21 7 7 686 28 14 49 28 21 14 14 35 35 49 14 21 7 7 28 49 28 28 7 119 14 21 21 14 21 7 147 84 7 7 21 7 35 119 49 7 14 7 21 14 7 7 28 28 7 21 14 14 7 98 7 7 14 21 28 35 119 119 14 7 14 14 7 35 14 7 21 7 7 21 63 7 7 91 14 56 7 14 42 14 21 14 28 28 7 77 7 21 14 7 21 7 7 7 35 35 7 91 14 14 7 7 14 21 21 14 28 70 14 217 70 322 14 28 84 7 7 7 21 49 7 14 7 14 14 91 56 70 7 14 21 28 7 70 28 56 21 777 35 42 28 14 21 7 28 70 42 7 14 7 63 14 70 140 28 28 7 14 7 7 77 28 14 7 7 7 147 7 7 35 14 21 7 14 14 14 63 14 70 14 35 28 7 35 35 105 49 56 21 154 14 119 7 42 14 154 70 28 14 14 14 63 7 98 7 14 56 28 35 7 7 112 42 7 28 7 7 14 14 7 14 21 35 154 56 49 7 112 35 14 14 7 28 7 7 7 7 7 14 7 70 35 126 42 56 21 49 21 49 35 35 14 21 14 28 21 28 14 7 7 7 21 14 35 14 35 133 7 7 21 21 35 42 14 7 42 98 777 7 42 49 7 21 28 28 7 7 35 14 49 98 35 7 7 14 7 14 35 21 14 14 56 28 7 119 7 7 28 21 7 28 28 7 84 147 14 14 35 7 7 7 7 175 70 7 7 7 70 21 21 21 119 7 7 42 7 14 7 7 7 28 105 21 14 63 7 7 7 63 7 112 49 42 70 56 105 35 21 21 28 7 14 7 28 49 14 28 28 7 91 56 21 7 28 7 28 35 224 14 28 49 14 7 7 14 63 7 7 7 35 28 21 14 28 7 21 14 7 91 7 7 28 14 35 7 28 7 140 7 21 56 7 35 56 28 42 7 21 7 21 14 35 14 42 21 35 77 14 287 7 98 7 35 35 7 21 14 7 7 28 42 7 14 56 7 28 91 49 14 14 14 14 7 7 42 35 42 7 98 21 14 14 7 35 7 56 14 56 7 35 7 42 28 7 98 21 7 7 14 49 21 70 14 7 21 21 7 7 28 14 28 161 7 28 14 77 7 7 7 7 7 21 14 21 21 294 84 245 7 28 56 28 7 49 14 14 35 70 7 63 7 7 7 49 49 7 42 7 28 28 21 14 7 21 77 7 42 14 14 21 112 14 434 21 56 14 7 21 35 28 14 14 7 35 7 7 14 35 7 231 14 14 7 42 28 196 7 70 7 14 70 7 7 7 112 7 7 14 63 14 49 42 196 21 14 14 7 42 77 14 7 21 21 7 35 7 28 70 7 7 14 14 21 28 56 28 175 119 7 21 133 7 14 21 7 133 21 49 21 14 35 7 105 14 84 42 14 105 14 28 7 21 63 35 14 35 84 49 28 105 7 7 7 14 14 7 84 21 56 21 21 84 14 28 14 14 14 7 63 49 7 7 91 70 28 56 21 98 63 7 28 21 14 14 21 35 28 49 21 7 42 7 63 35 14 28 14 35 28 7 35 14 14 21 42 21 182 35 7 7 21 14 21 7 28 14 35 14 315 35 7 175 28 28 7 189 14 14 7 28 14 7 7 7 91 7 14 14 154 7 49 21 14 119 14 28 28 7 42 28 63 7 7 14 42 21 21 21 7 14 21 42 7 21 119 14 21 14 7 21 7 14 84 21 14 14 21 14 7 84 28 21 21 28 7 147 14 98 7 56 14 14 49 112 7 7 133 7 63 35 224 28 7 147 35 21 14 35 14 112 21 308 91 14 49 42 133 21 7 7 35 14 14 7 28 7 7 35 91 182 21 14 7 56 28 14 14 49 7 28 35 140 49 42 21 35 7 28 7 21 7 7 14 35 224 14 14 35 7 14 7 91 7 28 7 14 7 42 7 98 14 7 7 28 7 7 42 7 28 14 35 7 7 91 28 28 56 42 84 21 42 266 56 280 35 7 7 28 7 7 133 7 77 14 14 21 35 7 14 126 91 42 28 56 14 7 7 14 35 7 28 14 147 14 7 7 42 14 126 77 14 42 7 49 35 28 14 56 77 70 7 14 42 7 133 14 7 7 63 7 112 70 7 7 7 35 63 84 49 70 28 133 35 56 14 49 7 133 56 147 42 42 28 21 49 21 70 7 91 21 35 7 7 7 7 21 14 14 7 35 7 56 7 7 21 42 7 14 14 56 238 14 7 7 7 7 7 133 7 21 7 7 7 28 21 7 28 14 21 28 35 14 154 7 7 140 21 168 49 56 28 14 7 42 7 14 651 14 119 14 28 14 7 42 14 42 21 63 56 77 35 7 21 35 14 21 63 7 21 21 35 28 14 14 175 112 14 91 14 77 28 14 14 21 21 56 154 28 14 560 28 21 91 7 21 7 28 35 7 7 42 21 21 35 7 28 21 77 7 7 7 35 14 42 14 56 343 14 35 7 70 147 7 147 21 21 42 14 7 21 28 7 7 21 14 91 35 49 84 7 77 140 84 35 28 161 14 21 7 7 7 70 7 7 28 7 42 21 14 7 14 14 63 14 7 7 161 7 7 7 28 14 77 21 63 7 7 14 28 301 133 63 21 14 21 35 56 294 175 21 168 70 7 28 14 35 14 14 70 14 7 56 7 7 49 35 175 7 21 7 28 7 42 21 7 84 14 56 7 7 7 42 42 7 84 21 7 7 84 112 7 7 14 119 49 84 98 7 35 28 14 14 35 7 35 77 28 14 21 28 28 7 42 7 28 7 14 7 7 14 63 14 7 14 35 7 28 63 7 105 35 7 21 252 7 14 21 35 14 7 21 21 14 7 21 14 14 84 7 63 35 7 7 63 70 56 70 7 7 7 35 7 84 126 49 28 21 21 14 7 7 7 7 84 35 14 56 35 7 14 21 42 14 28 7 14 28 28 7 7 21 14 7 7 7 182 14 14 7 105 28 7 84 35 42 21 56 14 14 7 126 7 70 49 70 28 56 14 21 70 21 21 35 28 259 14 35 7 35 161 7 14 49 42 70 21 280 7 35 14 7 7 7 21 35 21 154 35 84 7 42 21 70 70 21 14 14 14 7 70 7 63 49 14 21 14 70 7 21 42 14 126 49 21 112 7 28 385 14 14 7 14 21 14 49 35 7 63 91 7 7 14 14 14 35 7 203 21 56 35 7 98 14 49 7 7 7 7 28 7 147 7 7 70 14 70 14 21 42 7 91 56 70 21 14 21 35 14 14 294 28 14 7 119 7 7 7 7 14 21 7 7 42 56 28 35 35 14 70 21 7 7 7 7 21 35 91 28 7 21 7 7 49 21 42 7 77 7 7 7 21 84 7 7 7 7 35 777 14 7 28 7 28 56 28 56 7 98 28 7 7 14 70 21 35 84 7 21 49 7 91 7 21 21 63 77 42 7 56 56 7 7 63 42 7 49 7 14 7 14 49 14 7 7 105 7 49 7 7 7 14 35 35 7 21 7 14 7 7 21 7 7 21 28 98 28 28 28 21 21 7 14 49 21 14 49 49 7 7 14 35 42 14 49 7 14 7 63 14 21 49 7 7 7 7 77 28 14 21 7 28 7 42 56 7 294 7 105 35 14 21 28 14 7 28 14 14 35 14 28 49 14 7 77 14 63 42 70 203 7 42 63 119 14 7 7 21 7 21 7 70 105 21 14 7 28 28 63 42 168 35 35 21 28 21 14 14 56 28 168 56 147 42 7 21 14 28 70 532 7 7 203 14 42 70 7 63 7 42 35 7 189 84 91 35 21 14 14 49 28 7 28 77 7 14 21 7 14 7 7 21 21 7 28 7 21 7 14 14 7 154 21 7 7 7 7 21 35 7 7 35 35 70 7 28 14 28 7 21 14 42 56 7 56 14 49 28 378 7 21 14 35 14 7 28 98 105 14 14 35 105 28 28 7 7 14 35 7 56 7 175 91 7 14 21 35 7 35 14 21 21 42 7 7 14 35 35 21 7 14 14 14 14 7 21 7 21 21 21 63 28 49 14 35 14 7 7 91 7 14 245 14 112 56 7 63 28 63 77 14 70 63 7 35 63 42 21 49 35 28 35 14 21 7 7 14 28 63 105 28 49 7 14 7 14 7 21 14 539 77 7 42 14 14 63 7 56 49 35 42 7 28 112 91 7 28 14 14 21 21 35 7 14 154 28 56 7 56 77 35 14 7 105 49 14 28 21 56 217 35 21 7 28 49 56 154 21 126 7 84 70 70 98 49 63 98 14 14 14 7 63 14 35 7 63 21 7 7 42 133 7 14 14 84 7 105 49 28 21 7 21 7 28 14 42 7 7 14 49 28 133 28 7 21 7 252 14 28 7 35 35 21 7 49 7 28 7 35 77 77 7 7 77 63 49 105 126 7 70 21 189 14 7 14 63 28 14 14 21 7 49 91 7 7 154 42 7 14 119 7 77 98 21 49 56 7 21 14 14 14 7 14 14 14 7 28 7 28 7 7 14 14 49 14 203 7 49 7 7 7 7 7 28 21 49 56 28 7 21 35 7 14 7 7 21 28 7 35 35 7 7 287 7 7 14 42 56 49 21 42 7 21 7 7 7 56 42 98 7 42 14 21 14 14 14 63 21 28 35 28 14 7 7 7 14 7 105 21 42 49 35 7 21 14 56 56 28 21 203 7 14 7 7 7 7 49 7 7 21 21 14 147 14 35 21 7 14 35 56 21 35 7 777 28 84 14 315 14 7 14 7 21 35 7 7 14 161 266 35 7 7 7 7 7 14 14 21 21 21 42 21 21 7 7 14 7 140 210 21 28 21 14 21 7 14 35 21 14 77 70 21 7 56 7 42 14 14 7 14 63 14 7 70 21 42 21 21 126 231 14 7 14 42 84 14 49 35 28 14 14 14 7 77 91 14 7 7 7 21 28 7 84 28 21 7 14 14 21 7 252 14 28 7 21 14 7 14 322 7 7 98 77 14 21 7 245 7 56 35 14 21 63 14 35 7 7 14 21 28 14 42 105 119 63 42 14 42 63 14 21 7 14 7 14 91 14 21 21 7 14 14 7 14 14 70 14 777 14 91 28 14 7 7 7 56 14 7 126 7 7 21 7 14 7 7 14 14 14 7 7 14 7 28 14 28 21 7 7 7 35 7 119 14 7 7 21 35 42 91 14 21 7 7 7 28 91 28 14 7 7 84 21 7 56 28 21 14 35 77 70 7 7 21 28 49 7 105 28 49 7 14 7 119 42 7 14 28 7 7 28 7 14 7 49 7 182 14 35 77 35 14 7 175 7 21 112 7 14 7 21 14 7 7 7 7 14 63 28 7 42 14 7 14 14 140 63 7 21 35 35 7 14 14 21 84 63 35 14 7 14 14 21 28 182 35 14 14 35 21 77 21 7 35 7 42 42 14 7 7 28 28 21 77 7 14 14 7 7 7 14 7 56 84 7 28 21 49 28 56 7 7 21 28 77 14 21 84 7 7 21 7 35 21 49 7 28 21 7 126 21 28 14 7 7 28 7 7 7 49 112 224 7 35 7 49 42 28 49 28 14 42 21 14 91 28 14 14 7 168 49 21 7 49 7 119 98 14 7 7 49 7 91 28 42 7 7 14 14 56 126 14 14 56 7 14 56 14 7 63 42 28 14 56 49 21 77 7 21 21 14 35 7 7 112 238 35 35 7 70 28 70 7 196 7 14 14 21 42 7 14 21 7 21 35 28 35 7 105 105 7 35 7 7 7 42 7 14 28 21 28 21 14 7 7 70 7 7 56 7 21 21 7 77 7 7 7 28 35 21 35 7 7 7 133 56 7 7 21 21 21 14 196 28 42 7 7 35 7 77 7 7 7 14 14 28 42 7 14 21 14 133 21 7 35 14 7 21 21 70 7 7 133 7 7 21 21 21 14 7 42 35 28 21 70 14 7 70 14 42 7 21 196 7 21 21 7 14 49 21 7 21 7 7 21 35 7 7 21 42 42 14 49 21 7 28 28 14 7 7 14 28 7 14 14 84 14 42 14 49 14 7 28 42 7 7 21 35 21 28 7 21 77 28 14 42 7 42 21 28 14 7 14 42 14 7 7 7 21 91 7 7 21 35 21 42 63 21 7 7 84 14 84 126 14 42 14 7 7 91 14 7 14 28 7 28 14 28 7 56 42 7 133 7 56 14 21 658 7 14 42 77 63 35 63 7 7 7 7 7 21 28 49 28 7 42 21 7 14 7 21 7 7 49 28 7 28 14 14 70 14 7 7 7 21 21 28 42 63 7 77 7 28 7 28 105 21 21 7 7 63 21 35 42 98 14 7 42 42 21 105 49 14 14 28 21 7 7 21 322 21 35 7 35 14 14 7 133 7 105 119 35 7 98 35 14 21 98 14 21 49 21 14 14 7 7 7 91 7 7 56 49 42 63 7 56 70 84 7 14 133 7 49 7 7 28 14 119 21 21 7 14 7 7 42 14 7 14 7 14 42 21 7 14 7 7 7 7 28 14 7 42 427 14 49 7 7 21 14 336 56 56 21 7 7 42 28 28 7 7 7 161 14 28 14 28 7 56 7 42 14 112 21 7 77 7 28 42 28 7 14 14 21 14 49 14 7 70 70 14 98 14 14 161 7 14 14 7 35 77 49 21 7 21 7 21 56 112 14 21 91 28 14 21 7 28 21 14 14 56 35 7 147 28 14 7 7 7 7 28 7 56 7 63 42 35 21 7 7 343 49 175 7 70 21 42 14 21 56 28 14 35 35 14 98 7 63 7 35 21 70 119 63 126 7 14 28 49 49 14 35 70 28 7 21 91 7 7 21 7 7 7 7 21 7 14 21 28 21 14 14 28 168 77 7 7 7 7 7 406 21 7 14 112 7 77 119 49 49 14 7 21 14 21 14 14 49 154 21 7 259 7 112 7 21 7 7 28 63 21 14 7 14 28 28 14 7 70 21 42 42 7 21 35 7 7 28 119 56 7 7 21 70 14 7 21 14 28 42 28 7 7 7 49 21 14 70 42 14 133 7 49 14 49 28 35 14 14 7 35 105 14 105 112 21 7 7 35 21 7 7 14 21 21 7 147 7 28 14 63 14 21 35 7 21 7 7 49 14 14 7 42 7 14 14 21 7 21 21 14 56 42 7 266 21 35 7 154 14 14 14 154 21 35 14 7 63 14 14 42 7 7 7 14 14 49 7 7 14 21 14 539 21 189 28 7 7 364 7 28 21 7 21 84 21 161 7 14 7 7 7 14 63 35 140 7 56 98 154 21 35 7 28 7 14 7 7 35 7 35 14 7 7 35 14 7 28 28 42 35 7 14 42 119 7 14 35 28 14 7 28 7 35 7 7 70 35 7 21 42 28 7 7 14 7 14 28 42 28 14 119 70 28 14 14 21 49 14 14 21 196 63 63 7 7 49 7 14 7 252 42 7 14 14 14 63 42 14 7 7 42 7 14 7 63 133 126 7 126 14 35 126 7 14 7 28 42 119 70 7 21 7 42 14 63 14 28 14 49 49 21 7 7 7 35 14 28 14 7 7 105 63 14 21 49 7 28 35 7 7 7 35 14 28 42 35 238 7 7 7 7 42 14 133 7 42 7 21 42 14 14 140 35 21 28 7 21 7 14 42 21 84 7 7 112 21 91 14 7 42 28 7 21 21 7 28 63 77 42 14 7 98 7 7 42 7 28 42 21 28 28 14 98 105 42 21 28 21 21 21 7 7 28 7 49 7 7 14 7 14 14 35 14 42 84 35 42 56 21 14 35 112 126 21 28 133 7 7 7 7 140 14 7 21 7 7 161 7 14 7 7 98 147 21 35 14 7 7 14 21 28 14 7 56 63 7 21 21 7 7 63 28 7 91 7 21 21 63 28 28 21 7 35 7 14 147 14 42 7 7 14 56 28 7 7 14 21 14 42 42 35 21 7 7 56 7 35 14 210 14 14 7 7 42 28 7 56 49 7 21 7 161 35 21 28 7 84 14 14 7 49 21 7 7 98 7 21 7 105 14 14 350 7 56 28 21 7 91 14 84 21 7 7 21 14 35 14 14 14 49 14 28 42 14 14 21 147 35 14 98 7 21 42 56 42 7 42 7 14 7 14 7 42 14 7 7 7 7 7 42 21 7 14 14 14 63 35 49 21 42 35 7 28 28 70 14 14 119 21 21 14 7 7 77 14 175 7 28 28 105 28 21 14 7 14 21 98 28 35 21 42 7 14 14 28 7 7 28 147 21 14 21 14 7 70 14 21 7 14 84 28 14 7 7 56 7 42 7 77 77 63 35 7 28 56 21 63 7 14 70 14 14 7 14 35 21 14 21 14 14 14 14 28 28 21 56 28 7 21 7 35 42 7 14 7 14 14 7 189 7 21 21 14 42 42 28 42 7 14 7 14 77 7 14 7 35 14 7 42 7 21 7 77 42 21 42 7 7 217 35 7 189 7 21 21 21 14 28 28 14 14 49 14 28 14 7 7 56 7 56 7 35 14 7 28 84 7 14 28 7 35 21 126 42 14 7 42 35 56 7 14 21 105 49 35 7 35 7 21 21 14 42 28 28 42 14 7 21 21 217 7 7 140 14 21 7 28 14 84 42 14 21 14 14 7 7 28 14 28 14 21 70 35 49 7 14 28 49 28 7 14 70 49 42 42 35 7 7 28 14 21 35 21 28 7 63 35 70 7 21 49 42 7 14 14 7 77 28 35 7 7 7 14 21 161 7 14 686 14 14 28 7 21 70 14 14 21 35 63 7 21 147 21 14 21 21 14 14 28 21 28 28 7 49 70 7 7 14 7 7 7 7 14 49 35 7 42 14 14 35 14 7 7 7 7 14 14 7 21 7 119 14 35 35 14 7 28 7 14 63 7 7 21 35 7 28 7 21 7 56 98 7 7 42 7 14 434 7 7 126 7 7 28 14 21 7 56 308 63 140 14 14 56 14 21 35 35 7 7 7 21 28 7 105 7 21 161 35 7 21 154 56 35 168 7 49 49 133 35 14 7 14 7 14 21 14 28 7 14 7 14 7 28 14 7 7 7 7 7 28 14 28 14 7 7 21 7 7 56 7 14 119 7 21 21 14 28 14 14 84 287 21 84 84 49 7 7 7 98 7 63 42 35 21 35 21 14 28 7 14 7 28 14 7 21 28 84 28 7 119 7 7 7 189 28 21 7 7 7 28 7 7 28 91 7 49 14 7 63 7 49 42 49 28 14 252 14 21 35 7 21 112 21 133 49 28 14 35 21 14 7 7 14 21 42 7 7 210 84 7 35 42 28 21 126 14 28 203 91 35 63 42 28 21 98 49 7 14 21 273 14 7 28 266 7 14 28 21 56 7 14 35 14 7 7 14 105 126 42 7 98 35 49 63 14 49 21 28 21 7 21 35 28 7 42 35 266 21 21 14 14 7 7 49 49 14 7 77 7 56 7 21 112 7 21 35 63 28 42 42 28 7 714 14 21 35 35 7 63 7 63 49 14 7 42 35 7 7 35 7 77 49 7 21 7 21 7 70 7 56 7 42 14 49 14 7 28 28 42 21 35 49 28 14 14 7 7 7 14 7 84 21 14 21 14 42 14 42 7 7 63 28 14 21 49 14 14 7 7 49 21 28 35 14 14 49 49 7 7 7 7 35 7 42 7 7 7 28 14 119 14 133 14 7 14 91 343 14 14 42 7 7 7 14 7 7 7 14 28 7 35 14 14 14 7 63 91 168 14 49 7 7 91 56 56 7 21 147 21 14 21 7 7 35 7 14 91 35 7 35 42 7 406 7 7 7 21 7 49 133 14 21 98 14 77 14 7 35 21 63 70 112 21 14 35 21 21 7 70 42 21 56 7 49 14 21 91 7 21 7 7 42 14 28 14 70 84 21 14 21 7 56 63 14 91 7 35 7 14 7 21 42 56 70 28 28 35 14 14 105 49 154 7 7 7 7 42 21 7 70 14 7 49 14 7 14 21 63 35 14 7 56 28 14 14 7 7 7 49 70 56 28 14 42 7 161 14 14 14 35 7 7 42 175 28 14 7 42 119 7 21 63 140 105 14 28 49 91 35 14 7 14 7 7 42 21 7 28 70 56 14 35 7 14 35 7 7 14 7 21 21 350 147 7 7 14 49 14 7 28 14 21 14 28 42 28 35 14 21 7 7 7 14 21 14 7 21 42 91 7 7 161 7 35 7 21 14 21 98 105 42 7 21 42 14 7 35 14 91 28 14 7 105 14 35 42 28 21 84 14 21 7 21 21 7 7 7 21 7 7 70 28 21 35 777 35 196 7 7 98 21 35 7 42 21 7 21 7 14 21 42 21 28 21 77 28 7 84 7 35 21 7 70 7 7 7 245 231 21 175 7 28 7 7 7 161 7 56 7 14 14 49 7 7 28 14 7 7 7 77 7 7 56 7 42 14 21 7 56 14 21 413 84 7 28 91 7 7 7 7 469 21 7 7 294 14 7 14 63 56 42 14 7 14 84 35 35 42 7 7 56 14 7 14 7 35 14 21 42 7 7 70 21 7 140 7 7 49 7 7 49 112 49 42 77 77 21 14 21 7 42 14 35 7 7 7 7 7 14 14 196 133 42 49 14 7 91 133 14 7 7 98 77 21 112 7 77 7 7 21 35 7 28 14 77 7 7 7 21 7 63 35 21 133 77 7 28 133 14 7 7 42 63 14 112 21 91 7 14 21 7 28 84 21 28 28 42 21 28 7 63 133 21 77 77 7 7 21 98 28 35 21 7 28 35 14 35 21 7 7 154 63 21 14 7 112 28 42 84 42 14 7 28 7 14 14 49 7 35 182 14 7 112 63 7 14 378 7 42 7 7 7 7 7 7 21 7 49 7 14 14 7 49 28 14 7 21 7 91 21 7 56 21 7 7 7 7 7 14 14 147 14 14 7 28 14 28 7 14 7 7 21 7 21 14 14 98 91 14 14 21 56 35 7 35 7 42 28 7 21 21 35 7 7 7 14 49 14 28 7 56 7 35 28 28 7 56 63 7 70 49 7 14 21 7 14 98 42 14 21 21 7 35 7 98 42 147 7 42 42 42 28 7 7 14 7 7 21 14 28 35 28 14 49 7 119 56 21 28 35 112 14 119 35 35 14 7 35 7 14 14 28 84 28 35 21 7 77 14 42 14 84 14 98 21 7 21 21 35 14 7 14 28 35 7 28 7 14 21 7 7 28 35 7 7 28 7 28 21 14 175 7 49 56 210 28 56 28 7 35 14 21 7 161 35 7 7 28 63 21 14 14 7 14 63 91 63 7 7 14 84 28 175 7 14 7 203 7 7 7 28 35 35 7 35 7 7 14 28 35 49 91 35 28 7 14 28 21 14 7 35 35 14 28 21 14 7 21 14 434 7 63 14 7 7 14 70 14 14 28 14 7 7 7 35 21 28 14 7 21 42 35 7 77 14 7 35 7 14 14 63 28 7 56 21 21 14 42 7 42 7 21 42 7 21 7 14 14 7 14 70 679 7 21 7 133 28 7 21 21 7 35 49 7 21 7 7 35 35 7 7 63 7 28 35 7 28 28 119 7 7 21 7 133 42 14 28 7 56 14 14 7 7 14 7 14 28 21 7 35 14 14 7 14 105 28 91 14 35 49 28 21 14 238 7 28 7 56 21 7 140 14 21 70 49 7 42 70 7 28 42 7 14 7 49 21 7 14 21 7 14 28 7 7 35 7 14 21 21 28 14 35 14 112 49 42 371 28 14 7 42 14 28 14 7 21 21 7 14 28 56 175 21 42 7 21 7 42 14 21 56 7 49 14 28 7 56 7 21 35 28 98 7 98 21 56 7 49 28 245 77 7 14 14 7 7 21 147 7 14 49 21 49 42 7 28 14 7 7 28 7 7 7 7 35 28 21 56 7 35 77 105 14 14 7 7 7 7 14 49 35 14 112 7 56 28 14 21 42 21 7 7 7 14 49 7 7 308 7 14 126 63 28 14 7 7 28 14 42 7 14 7 21 63 7 56 28 63 77 7 14 7 14 7 14 14 7 14 7 91 70 28 42 7 21 28 63 175 14 35 7 7 49 7 56 14 14 84 56 21 35 21 21 42 21 14 7 21 35 42 35 7 14 63 14 14 35 7 28 21 28 49 7 14 49 14 14 35 7 7 42 7 7 35 14 21 42 42 14 35 42 14 21 7 56 14 49 7 7 7 28 7 28 21 7 7 7 14 35 7 28 7 28 77 21 14 14 49 21 7 14 21 28 7 14 7 14 126 28 168 7 21 133 49 7 7 14 49 14 7 21 7 14 28 7 35 28 14 14 7 42 7 7 7 21 28 21 28 28 35 7 49 56 126 231 21 7 14 14 56 7 21 77 49 21 7 7 14 28 14 21 189 7 154 14 28 70 7 14 14 21 7 21 21 7 21 7 14 119 35 14 119 63 14 42 28 14 84 77 7 28 42 49 7 14 14 7 182 35 238 42 42 7 133 98 14 28 63 14 42 14 7 7 7 7 42 7 14 14 7 84 21 14 7 49 175 7 56 14 7 14 77 21 91 7 28 7 7 70 28 7 98 98 42 21 14 21 7 70 98 35 14 7 7 35 35 112 35 7 77 7 84 63 7 70 49 14 7 14 42 49 7 28 126 56 56 28 147 7 14 7 21 21 14 161 14 28 84 42 413 21 14 7 7 98 7 84 42 21 98 42 7 28 7 35 21 49 21 14 21 21 35 14 14 7 7 14 7 28 7 119 28 21 14 28 28 7 21 14 21 14 14 7 14 56 49 21 7 49 7 7 49 7 14 77 28 56 35 28 7 28 21 28 28 7 14 14 42 42 21 21 42 63 42 28 133 77 35 56 63 7 21 21 28 35 35 14 28 49 7 63 14 49 70 35 28 77 14 70 21 112 7 14 21 7 14 7 70 14 35 7 7 7 7 28 21 7 49 21 21 7 35 77 553 245 21 14 49 49 21 14 14 35 7 28 98 14 7 56 21 28 175 7 14 14 7 21 7 14 7 21 21 14 14 91 35 7 21 7 14 21 77 35 14 7 35 35 7 21 49 7 21 154 7 63 35 7 7 21 14 7 42 42 28 28 56 63 147 14 7 14 175 98 42 28 35 21 7 21 21 91 7 14 63 7 21 14 7 14 49 7 21 14 35 7 77 28 7 21 7 112 7 14 21 7 14 77 56 7 28 7 7 28 14 21 14 7 91 28 28 147 7 14 28 49 21 7 35 14 7 63 28 7 224 7 42 21 7 7 42 119 84 7 7 7 7 21 7 7 28 224 7 63 42 21 35 7 7 28 7 28 14 14 14 14 35 14 49 21 7 28 28 7 35 7 35 28 14 21 7 7 378 7 63 7 21 77 14 581 14 7 35 28 7 35 14 14 28 126 7 7 35 7 21 301 63 7 147 63 42 7 7 21 21 7 14 119 14 21 7 7 14 7 21 7 28 140 28 49 21 14 7 7 42 28 21 21 14 7 14 28 7 28 42 7 322 84 42 7 28 14 42 7 7 14 7 28 14 14 7 14 49 168 14 7 7 28 189 7 49 21 14 63 14 7 14 28 14 14 7 14 14 35 413 119 63 21 35 7 56 7 7 28 56 49 21 35 91 28 7 84 21 7 35 7 7 49 42 14 7 14 14 14 49 70 35 7 7 49 21 49 7 14 7 7 21 7 70 21 84 21 28 7 7 7 7 119 63 56 28 14 49 14 7 14 21 7 7 21 7 14 28 14 70 21 7 49 21 7 7 42 7 14 14 7 63 28 7 14 21 14 7 7 98 7 14 14 42 28 7 98 7 21 7 91 21 14 77 35 14 21 14 7 21 7 35 14 42 126 21 28 14 7 77 56 35 42 7 7 14 21 7 14 14 42 21 56 7 7 84 56 14 7 14 14 98 63 28 7 7 14 7 56 14 42 7 14 77 14 42 21 70 42 7 21 21 7 21 105 63 7 7 21 98 7 7 21 63 63 14 14 42 14 7 119 14 7 42 14 21 21 42 777 7 56 7 42 14 28 91 14 7 98 42 14 7 7 14 7 7 7 7 14 28 14 35 14 14 21 28 7 21 7 42 434 14 28 105 14 28 84 7 42 7 42 7 21 21 7 7 56 14 98 7 7 21 7 21 42 35 14 42 7 182 7 70 35 35 42 7 14 42 56 21 28 42 28 14 7 7 14 49 35 224 91 21 42 21 35 7 14 7 21 21 49 14 7 14 7 21 7 70 14 14 7 7 56 21 14 91 35 7 35 7 35 14 49 21 21 84 42 28 49 7 84 7 7 7 21 35 21 21 14 7 84 21 49 7 42 42 7 7 35 84 28 14 21 14 56 49 7 7 14 28 28 7 14 7 21 14 7 21 56 42 14 49 98 7 7 7 105 175 28 21 7 14 7 14 21 175 70 42 21 28 21 126 28 7 7 77 7 21 21 7 14 35 7 28 14 21 14 28 14 84 7 14 21 224 35 28 350 14 42 21 7 14 14 35 7 7 21 28 7 238 112 21 133 7 28 28 35 7 28 21 28 154 63 336 7 161 14 7 35 42 7 14 7 21 56 70 28 35 28 63 70 7 21 7 14 28 7 126 21 14 140 35 7 56 7 28 7 7 21 63 217 42 14 21 14 7 7 14 35 7 14 28 7 7 21 7 49 7 35 21 21 21 7 28 28 210 35 28 133 14 7 14 7 21 49 7 14 21 7 14 259 56 7 63 231 28 7 546 21 7 56 14 7 28 21 28 77 7 21 35 14 42 14 35 7 7 14 77 14 14 35 105 7 14 21 28 7 21 7 35 21 7 14 21 35 7 14 7 35 266 7 21 7 35 28 84 133 7 7 35 49 14 14 70 56 21 14 21 7 14 14 28 28 7 28 42 7 14 28 49 21 42 84 7 63 119 21 77 14 14 7 21 28 42 14 7 49 28 7 7 70 49 56 14 28 14 7 70 7 14 7 28 42 154 42 21 35 7 112 14 56 7 7 42 154 7 336 35 98 21 28 14 56 7 21 56 7 77 63 7 7 28 14 112 21 14 28 21 35 14 14 35 35 7 7 7 35 21 7 14 7 49 91 14 42 133 14 7 7 14 35 14 21 21 7 28 14 7 7 7 7 42 7 14 77 7 35 21 14 49 7 21 21 28 259 7 14 7 7 28 14 42 21 14 35 133 63 28 28 7 28 28 28 14 14 14 35 14 7 7 84 91 21 14 14 42 42 63 28 28 7 42 28 14 28 7 280 7 28 28 91 7 7 28 21 7 7 28 252 14 21 28 7 14 7 28 49 63 21 7 14 35 56 35 28 21 35 49 14 63 14 7 14 133 49 7 21 28 7 7 14 7 7 28 21 42 7 14 105 28 84 98 7 35 42 14 28 21 77 7 28 7 21 112 21 7 35 14 7 70 28 21 7 7 21 7 28 35 7 7 14 14 63 7 672 21 63 14 7 273 21 35 42 84 14 42 77 14 7 14 245 7 154 35 14 7 7 35 35 7 7 77 7];
% countconn = zeros(1,maxconn-minconn+1);
% for i = 1:10000
%     countconn(actualconn(i)/7) = countconn(actualconn(i)/7) + 1;
% end
% 
% figure;
% axes('FontSize',24);
% plot(minconn:maxconn,countconn./10000);grid on;title('Connections vs Probability, Simulation');xlabel('Number of connections');ylabel('Probability of connections');