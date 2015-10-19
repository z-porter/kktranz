function rechi=kkrebook2z(omega,imchi,alpha, omegalow, omegahigh)
%Calculates real(chi) based on KK transform of input
%
%The program inputs are:
%the vector of the frequency (or energy) components,
%the vector of the imaginary part of the susceptibility under examination,
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
%If imchi is the imaginary part of a linear susceptibility,
%alpha must be 0.
%
%If imchi is the imaginary part of the nth
%harmonic generation susceptibility, alpha=0,1,..2n.
%
%If imchi is the imaginary part of a pump and probe
%susceptibility, alpha=0 or 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These files accompany the book
%"Kramers-Kronig Relations in Optical Materials Research"
%by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M.
%Springer, Heidelberg, 2005
%where the theory and applications are fully developed.
%The output is the estimate of the real part as obtained
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
end;
if size(imchi,1)>size(imchi,2);
    imchi=imchi';
end;
%Here the program rearranges the two vectors so that,
%whichever their initial shape, they become row vectors.
%}

g=size(omega,2);
%Size of the vectors.%
rechi=zeros(size(imchi));
%The output is initialized.
a=rechi; b=a;
%Two vectors for intermediate calculations are initialized
deltaomega= diff(omega(2:end));
%Here we compute the frequency (or energy) intervals

b(1)=sum(imchi(2:g).*omega(2:g).^(2*alpha+1)./(omega(2:g).^2-omega(1)^2));
rechi(1)=2/pi*deltaomega(1)*b(1)*omega(1)^(-2*alpha);
%First element of the output: the principal part integration
%is computed by excluding the first element of the input

a(g)=sum(imchi(1:g-1).*omega(1:g-1).^(2*alpha+1)/(omega(1:g-1).^2-omega(g)^2));
rechi(g)=2/pi*deltaomega(end)*a(g)*omega(g)^(-2*alpha);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input

%Specify subrange of omega values to report
if nargin <= 3
	omegalow = omega(1); 
	omegahigh= omega(end);
end
[~,idxL,~]  = find(omega > omegalow);
[~,idxH,~] = find(omega < omegahigh);
idxSubrange = intersect(idxL,idxH);

for j=idxSubrange
    a(j)=sum(imchi(1:j-1).*omega(1:j-1).^(2*alpha+1)./(omega(1:j-1).^2-omega(j)^2));
    b(j)=sum(imchi(j+1:g).*omega(j+1:g).^(2*alpha+1)./(omega(j+1:g).^2-omega(j)^2));
end

rechi(2:g-1) = (2/pi).*deltaomega.*(a(2:g-1)+b(2:g-1)).*omega(2:g-1).^(-2*alpha);
%Last element of the output: the principal part integration
%is computed by excluding the last element of