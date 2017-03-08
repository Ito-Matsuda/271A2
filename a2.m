% Written by Jose Manuel Matsuda 10150406
% I verify that this is of my own work

%https://www.mathworks.com/matlabcentral/fileexchange/24322-principal-component-analysis--pca--in-matlab/content/pca.m
function [] = a2()
%For me, make sure you are in Z: Desktop so you can load these files, or
%just copy the dat files to wherever the default to this thing is
load('z1.dat'); % This is 1a, already in matrix
load('z2.dat');
[evalvec,meanvec,evecmat] =  pcaprelim(z1); % Perform PCA and put em into variables
% 
% new_data = evalvec?
% step 2 is pcaaprox
%determine how many principal components are needed to capture most of
%variaton from mean signal (try with 5?)
% [approxcomp,approxvec]=pcaapprox(new_data,approxnum,meanvec,evecmat);
% must obtain co-ords of data points in direction of eigen vectors

% Using data signal 1, approxnum 2, and the data from pcaprelim
[approxcomp, approxvec] = pcaapprox(z2(:,21),2,meanvec,evecmat);
hold all % Used to plot on same chart
plot(approxvec,'--'); % Our approximation represented by dotted line reconstructed
plot(meanvec);
plot(z2(:,21));
% TODO, calculate the reconstruction errors
% Taking code from the Professor's notes
% http://research.cs.queensu.ca/~cisc271/pdf/class18.pdf

function [evalvec, meanvec, evecmat] = pcaprelim(Z)
% FUNCTION [EVALVEC, MEANVEC, EVECMAT] = PCAPRELIM(Z_DATA)
% performs the preliminary Principal Components Analysis
% (PCA) of A, a matrix in which the data are
% represented as columns. 
% PCAPRELIM returns:
% EVALVEC - eigenvalues of the PCA, sorted in decreasing order
%   holds most covariance
% MEANVEC - the mean vector of the initial data
%   
% EVECMAT - the eigenvectors of the PCA, as column vectors
%  ^ is 1e 
% Find the mean vector and form it into a matrix
[m,n] = size(Z); % M is the row size, n is the column size 
meanvec = mean(Z,2); %1b mean signal
% Takes the mean of each of the rows and puts it into a vector 
M = meanvec*ones(1,n); %1b mean matrix 
% Makes the vector a 101*30 matrix (same size as og) 

% Find the difference matrix and the covariance matrix
D = Z - M; % Difference matrix, Original - step b

% The co-variance matrix - step c
C=D*D'/(m-1); 

% Find the eigenvectors as a matrix and
% the eigenvalues as a vector step d
[Vvecs, Vvals] = eig(C); 
% ^ is Diagonal matrix Vvals of eigenvalues
% and matrix Vvecs whose columns are corresponding right eigenvec
% so that C*Vvecs = Vvecs*Vvals
% Calculating the co-effs of the principal components and their respective
% eigenfunctions from the co-variance mat
% Vvecs has the co-effs for principal components

% Diag elements of D store the variance of respective principal components
Vdiag = diag(Vvals); 
% Creates a single column containing the eigenvalues on the diagonal from Vvals

% Sort eigenvals from largest to smallest
[temp, Vval_index] = sort(Vdiag, 1, 'descend'); 

% Sort the original decomposition into the return variables
evalvec = Vdiag(Vval_index); %Eigenvalues of the PCA (aka covariance)
% that which capture the most variance 1-e

evecmat = Vvecs(:, Vval_index); % The Vval_index column of Vvecs
% Creates a matrix based off Vvecs using the columns values in Vval_index
% Eigenvecs of PCA as col vects (aka covariance)  1-e

% SELECTING A VALUE K BASED ON THE EIGENVALUES? start off with 5? lmao


%Again from the class notes

function [approxcomp,approxvec] = pcaapprox(new_data,approxnum,meanvec,evecmat)
% approximates new data based on a Principal Components Analysis
% (PCA) of initial data. Inputs are:
% NEW_DATA - a column vector to be approximated -> the column whose data
% will be used to plot it? use any of the columns of z1 say z1(:,3) the
% third column

% APPROXNUM - a scalar giving the order of the approximation 
    %should this just be 2 or 3?
% MEANVEC - the PCA mean vector (from PCAPRELIM)
    % use this to plot against? so you can see the difference? ie) the OG
% EVECMAT - the eigenvectors of the PCA(from PCAPRELIM)
%
% Return values are:
% APPROXCOMP - the components as a row vector of scalars
% APPROXVEC - the approximation of the new data as a vector
% Set up the return values

%Difference vector 
diffvec = new_data - meanvec;

approxcomp = zeros(approxnum, 1); 
% zeroes matrix with approxnum rows and 1 column
approxvec = meanvec; % The PCA mean vector from pcaprelim

% Loop through the eigenvectors, finding the components
% and building the approximation
for i=1:approxnum % From 1 to number of rows
evec = evecmat(:,i);
beta = dot(diffvec, evec);
approxcomp(i,1) = beta; % 2c
approxvec = approxvec + beta*evec; % 2d
end

%Plotting the mean signals, principal components, and reconstruction errors
%Dont plot all 60 signals, 