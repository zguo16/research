clear 
close all
clc


%% Generate the Shepp-Logan function

P = phantom( 'Modified Shepp-Logan',100 );% Generate the Shepp-Logan image
%P = normrnd(0,1,100,100);
imshow( P )


%% Set the parameter for numerical evaluation of FB expansion coefficients

c = 0.5;% bandlimit c

R = floor(size(P, 1)/2);% image lies in a disc of radius R

n_r = ceil(4*c*R);% number of points chosen in [0,c]

num_pool = 10;%?I do not understand that ,just copy that


%% The function computes computes the bessel functions
% The result after precomp_fb is:
%          basis.Phi_ns: bessel radial function J(R_{kq}\xi/c) for \xi in [0, c]
%          basis.ang_freqs: the associated angular frequencies (k in the paper)
%          basis.n_theta: the number of samples on concentric rings for polar Fourier transformation.
%          sample_points.r: positions in interval [0, c]
%          sample_points.w: weights

% ?I do not know the mean of weight here
[ basis, sample_points ] = precomp_fb( n_r, R, c );


%% Compute the Fourier-Bessel expansion coefficients

%Description of jobscript_FFBsPCA£º

% Computes Fourier-Bessel expansion coefficients 'coeff', and filtered steerable PCA expansion coefficients 'sPCA_coeff'.
%Input:
%	data: image dataset.
%	c: band limit
%	R: compact support radius
%	noise_variance: estimated noise variance
%Output
%	timing: toc_FFBsPCA: the whole time spent on computing fast steerable PCA. toc_sPCA: time for steerable PCA. toc_FBcoeff: time for computing Fourier-Bessel expansion coefficients
%	coeff: Fourier-Bessel expansion coefficients
%	mean_coeff: mean of the Fourier-Bessel expansion coefficients
%	sPCA_coeff: steerable PCA expansion coefficients
%	U: eigenvectors of C^{(k)}'s
%	D: eigenvalues of C^{(k)}'s. 

noise_variance = 0;% Additive Gaussian Noise

[ timing, coeff, mean_coeff, sPCA_coeff, U, D ] = jobscript_FFBsPCA(P, R, noise_variance, basis, sample_points, num_pool);


%% Generate the nonuniform distribution

around_1_degree = 0;% get the degree which is closest to 1, since we want to choose theta at 0,1,2,...,359

num_of_interval = 0;% number of interval we need to get around_1_degree.

delta_theta = 2*pi/basis.n_theta; % smallest interval we can choose.

while around_1_degree < 1./360*2*pi
    
    around_1_degree = around_1_degree + delta_theta;
    
    num_of_interval = num_of_interval + 1;
    
end

num_of_possible_theta = ceil(2*pi/around_1_degree);% the number of possile theta we can sanple.

pmf_theta = round(rand(num_of_possible_theta,1)*basis.n_theta);

pmf_theta = pmf_theta./sum(pmf_theta,1);%get the nonuniform pmf

alphabet = around_1_degree*(0:1:num_of_possible_theta-1)';%get the alphabet of theta

num_of_sample = 100;% the number of samples we want

rand_theta = randsrc(num_of_sample,1,[alphabet';pmf_theta']);%generate the random theta


%% Generating the sample theta

alpha = 30/360.*2*pi;% we want to cover the range of -30 to 30

num_of_delta_alpha = ceil(alpha./delta_theta);% number of delta_alpha we can have between 0 and 30

sample_degree = zeros(num_of_sample,2*num_of_delta_alpha+1);% record the sample. 

for i = -num_of_delta_alpha:1:num_of_delta_alpha
    
    sample_degree(:,i+num_of_delta_alpha+1) = rand_theta + delta_theta.*i.*ones(num_of_sample,1);
    
end

sample_degree = mod(sample_degree,2*pi);% mod 2*pi


%% Generate the projections
% do I need to include the normalize factor below?

k_max = size(coeff,1) - 1;%generate the k_max

P_theta_xi = zeros(num_of_sample,2*num_of_delta_alpha+1,n_r);% generate the projection

temp = zeros(1,1,n_r);

for i=1:1:num_of_sample
    
    for j=1:1:2*num_of_delta_alpha+1
        
        for k = -k_max:1:k_max
            
            if k < 0
                
                temp(1,1,:) = ((-1).^abs(k).*basis.Phi_ns{-k+1}.*exp(1i*k*sample_degree(i,j)))*conj(coeff{-k+1});
                
                P_theta_xi(i,j,:) = P_theta_xi(i,j,:) + temp;
                
                
            else
                
                temp(1,1,:) = (basis.Phi_ns{k+1}.*exp(1i*k*sample_degree(i,j)))*coeff{k+1};
                
                P_theta_xi(i,j,:) = P_theta_xi(i,j,:) + temp;
                
                
            end 
            
        end  
        
    end
    
end


%% Compute the covariance matrix

size_of_C = n_r*(2*num_of_delta_alpha + 1);

C_rho1_rho2_alpha1_alpha_2 = zeros( size_of_C, size_of_C );% Covariance matrix

 temp1 = zeros(num_of_sample,n_r);
 
 temp2 = zeros(num_of_sample,n_r);
 
for i = 1:1:(2*num_of_delta_alpha + 1)
    
    temp1(:,:) = P_theta_xi(:,i,:);
    
    for j = 1:1:(2*num_of_delta_alpha + 1)
        
       temp2(:,:) = P_theta_xi(:,j,:);
       
       C_rho1_rho2_alpha1_alpha_2(n_r*(i-1)+1:n_r*i,n_r*(j-1)+1:n_r*j) = temp1'*conj(temp2);
        
    end
    
end

 C_rho1_rho2_alpha1_alpha_2 = 1./num_of_sample.* C_rho1_rho2_alpha1_alpha_2;

 
 %% Generating the Psi matrix
 
 num_of_akq = size(basis.rad_freqs,1);
 
 Psi = zeros(size_of_C, 2.*(num_of_akq)-size(coeff{1,1},1));
 
 for i=-num_of_delta_alpha:1:num_of_delta_alpha
     
     temp3 = 0;
     
     for j = -k_max:1:k_max
         
         if j<0
             Psi((i+num_of_delta_alpha)*n_r+1:(i+num_of_delta_alpha+1)*n_r,temp3+1:temp3+size(coeff{abs(j)+1,1},1)) = ...
                 (-1).^abs(j).*basis.Phi_ns{abs(j)+1,1}.*exp(1i*j*i*delta_theta.*i);
         else
             Psi((i+num_of_delta_alpha)*n_r+1:(i+num_of_delta_alpha+1)*n_r,temp3+1:temp3+size(coeff{j+1,1},1)) = ...
                 basis.Phi_ns{j+1,1}.*exp(1i*j*i*delta_theta.*i);
             
         end
         
         temp3 = temp3 + size(coeff{abs(j)+1,1},1);
        
     end
     
 end
 
 
 %% Generating the P(\hat(C)) matrix
 
 %?? not sure about that. ask questions.
 
 FFT_pmf = fft(pmf_theta); % DFT of pmf
 
 C_P = zeros(2.*(num_of_akq)-size(coeff{1,1},1),2.*(num_of_akq)-size(coeff{1,1},1));
 
 temp4 = 0;
 
 for i = -k_max:1:k_max
     
     temp5 = 0;
     
     for j = -k_max:1:k_max
         if abs(i-j)<=size(FFT_pmf,1)-1
             C_P(temp4+1:temp4+size(coeff{abs(i)+1,1},1),temp5+1:temp5+size(coeff{abs(j)+1,1},1)) = ...
              FFT_pmf(abs(i-j)+1,1);
         else 
             C_P(temp4+1:temp4+size(coeff{abs(i)+1,1},1),temp5+1:temp5+size(coeff{abs(j)+1,1},1)) = 0;
         end
         
         temp5 = temp5 + size(coeff{abs(j)+1,1},1);
         
     end
     
     temp4 = temp4 + size(coeff{abs(i)+1,1},1);
     
 end
 
 
 
 
 







