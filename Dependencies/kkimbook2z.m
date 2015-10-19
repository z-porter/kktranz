function imchi=kkimbook2z(omega,rechi,alpha, omegalow, omegahigh)
%Calculates imag(chi) based on KK transform of input
%
%The program inputs are:
%the vector of the frequency (or energy) components,
%the vector of the real part of the susceptibility under examination,
%and the value of the moment considered.
%
%omegalow and omegaigh are optional and may be used to limit the omega
%region that is transformed (to save on computing time),
%but the output will still contain the entire input E range, with 0 for
%values below omegalow and above omegahigh.
%
%The two vectors must have the same length
%and the frequency vector omega must be equispaced.
%If not, apply MATLAB functions such as interp.
%
%If rechi is the real part of a linear susceptibility,
%alpha must be 0.
%
%If rechi is the real part of the nth
%harmonic generation susceptibility, alpha=0,1,..2n.
%
%If rechi is the real part of a pump and probe
%susceptibility, alpha=0 or 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This files accompanies the book
%"Kramers-Kronig Relations in Optical Materials Research"
%by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M.
%Springer, Heidelberg, 2005
%where the theory and applications are fully developed.
%The output is the estimate of the imaginary part as obtained
%with K-K relations.
%
%This software is distributed under the GNU licence agreement
%by Valerio Lucarini
%email: lucarini@alum.mit.edu
%University of Camerino
%Department of Mathematics and Computer Science
%Camerino, Italy
%
%vectorized/optimized by Zach Porter 7/30/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
if size(omega,1)>size(omega,2);
    omega=omega';
end; if size(rechi,1)>size(rechi,2);
    rechi=rechi';
end;
%Here the program rearranges the two vectors so that,
%whichever their initial shape, they become row vectors.
%}

%warning('step by step');
g=size(omega,2);
%Size of the vectors
imchi=zeros(size(rechi));
%The output is initialized.
a=imchi; b=a;
%Two vectors for intermediate calculations are initialized
deltaomega= diff(omega(2:end));
%Here we compute the frequency (or energy) intervals

b(1)=sum(rechi(2:g).*omega(2:g).^(2*alpha)./(omega(2:g).^2-omega(1)^2));
if isnan(b(1)); b(1) = 0;end;
imchi(1)=-2/pi*deltaomega(1)*b(1)*omega(1)^(1-2*alpha);
%First element of the output: the principal part integration
%is computed by excluding the first element of the input

a(g)=sum(rechi(1:g-1).*omega(1:g-1).^(2*alpha)/(omega(1:g-1).^2-omega(g)^2));
if isnan(a(g)); a(g) = 0;end;
imchi(g)=-2/pi*deltaomega(end)*a(g)*omega(g)^(1-2*alpha);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input.

if nargin <= 3
	omegalow = omega(1); 
	omegahigh= omega(end);
end
[~,idxL,~]  = find(omega > omegalow);
[~,idxH,~] = find(omega < omegahigh);
idxSubrange = intersect(idxL,idxH);

for j=idxSubrange
    a(j)=sum(rechi(1:j-1).*omega(1:j-1).^(2*alpha)./(omega(1:j-1).^2-omega(j)^2));
    b(j)=sum(rechi(j+1:g).*omega(j+1:g).^(2*alpha)./(omega(j+1:g).^2-omega(j)^2));
    %vectorizing is actually slower for small data sets
end

imchi(2:g-1) = -(2/pi).*deltaomega.*(a(2:g-1)+b(2:g-1)).*omega(2:g-1).^(1-2*alpha);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input