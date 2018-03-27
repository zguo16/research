function [ z ] = TomoRandRecov(coeff, basis, pmf_theta, C,num_of_delta_alpha,delta_theta,n_r)
%Description, solving quadratic equations
%Input:
%    coeff: expansion coefficients
%    basis: Fourier-Bessel Basis
%    pmf_theta: probability distribution of the views in theta
%    C covariance matrix constructed by projections!!
%Output:
%    z: recovered expansion coefficients. 
%Check if we can solve min \| \Psi ((aa').*C_P)) \Psi^* - C \|_F^2

%The following part is largely borrowed from shep_logan.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_max = size(coeff,1) - 1;%generate the k_max
%% Generating the Psi matrix (modified from shepp_logan.m)
 num_of_akq = size(basis.rad_freqs,1);
 %Psi = zeros(size_of_C, 2.*(num_of_akq)-size(coeff{1,1},1));
 
 for i=-num_of_delta_alpha:1:num_of_delta_alpha
     temp3 = 0;
     for j = -k_max:1:k_max
         if j<0
             %Psi((i+num_of_delta_alpha)*n_r+1:(i+num_of_delta_alpha+1)*n_r,temp3+1:temp3+size(coeff{abs(j)+1,1},1)) = ...
                 %(-1).^abs(j).*basis.Phi_ns{abs(j)+1,1}.*exp(1i*j*i*delta_theta);
              Psi((i+num_of_delta_alpha)*n_r+1:(i+num_of_delta_alpha+1)*n_r,temp3+1:temp3+size(coeff{abs(j)+1,1},1)) = ...
                 (-1).^abs(j).*fliplr(basis.Phi_ns{abs(j)+1,1}.*exp(1i*j*i*delta_theta)); %modified the order of the basis such that a_minus is easier to construct.
         else
             Psi((i+num_of_delta_alpha)*n_r+1:(i+num_of_delta_alpha+1)*n_r,temp3+1:temp3+size(coeff{j+1,1},1)) = ...
                 basis.Phi_ns{j+1,1}.*exp(1i*j*i*delta_theta); 
         end
         temp3 = temp3 + size(coeff{abs(j)+1,1},1);
     end 
 end

%% Generating the P(\hat(C)) matrix (directly copied from the shepp_logan.m)
 FFT_pmf = fft(pmf_theta); % DFT of pmf
 C_P = zeros(2.*(num_of_akq)-size(coeff{1,1},1),2.*(num_of_akq)-size(coeff{1,1},1));
 temp4 = 0;
 
 for i = -k_max:1:k_max    
     temp5 = 0;
     for j = -k_max:1:k_max
         if abs(i-j)<=size(FFT_pmf,1)-1
             if i-j<0
                 C_P(temp4+1:temp4+size(coeff{abs(i)+1,1},1),temp5+1:temp5+size(coeff{abs(j)+1,1},1)) = ...
                   FFT_pmf(abs(i-j)+1,1);
             else
                 C_P(temp4+1:temp4+size(coeff{abs(i)+1,1},1),temp5+1:temp5+size(coeff{abs(j)+1,1},1)) = ...
                   conj(FFT_pmf(i-j+1,1));
             end
         else 
             C_P(temp4+1:temp4+size(coeff{abs(i)+1,1},1),temp5+1:temp5+size(coeff{abs(j)+1,1},1)) = 0;
         end
         temp5 = temp5 + size(coeff{abs(j)+1,1},1); 
     end    
     temp4 = temp4 + size(coeff{abs(i)+1,1},1);    
 end
 
%% Get a (modified from shepp_logan.m)
%a1 = cell2mat(coeff);
%a_minus = cell2mat(coeff(k_max+1,1));
%a_minus = conj(a_minus);
%for i =k_max:-1:2
%    temp10 = cell2mat(coeff(i,1));
%    a_minus=[a_minus; conj(temp10)];
%end
%a=[a_minus;a1];

% Modify the structure for Psi, such that $a$ is easier to construct
% from the positive part
a1 = cell2mat(coeff);
p0 = size(coeff{1}, 1); %number of zero angular frequency expansion coefficients
a_minus = conj(flipud(a1(p0+1:end)));
a = [ a_minus; a1 ];
%C = Psi*((a*(a)').*C_P)*(Psi)';%compute the covariance matrix C 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve the optimization problem using minf_lbfgs or minf_lbfgsdl
z0 = randn(size(a1)) + sqrt(-1)*randn(size(a1)); %Because of symmetry, we just need ot initialize the positive coefficients. random initialization, need to think about a better way to initialize
F = @(z)objfun(z, C, Psi, C_P, p0);
options.TolFun = 1e-20;
options.TolX = 1e-20;
options.MaxIter = 50000;
%z = minf_lbfgs(F, [],z0,options);
z = minf_lbfgsdl(F, [],z0,options);
%z = [mean_sig*d; z];

end

%Define the cost function $\|\Psi (aa').*C_P \Psi^* - C \|_F^2
function [ cost ] = objfun( z, C, Psi, C_P, p0 )

z_minus = conj(flipud(z(p0+1:end))); %negative angular frequency coefficients
z = [ z_minus; z ]; %add the negative angular frequency coefficients.
Ca = Psi*(z*(z)'.*C_P)*(Psi)';
F = Ca - C;
cost = norm(F, 'fro')^2;

end
