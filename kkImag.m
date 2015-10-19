function [imf, varargout] = kkImag(ref, E_given, ElementName,...
    E_ROI_bounds,  ifplot, subtractZ)
%function for calculating f' from an input f'' using KK Transforms
%I call (f'=real[f]-f_0) ref and I call (f''=imag[f]) imf
%
%This will look choppy if ref does not match with table values at the ends,
%because of the KK Transform, since I do not shift ref in any way.
%
%Energy units are in eV
%
%INPUTS
%imf, an array of f'' data
%
%E_given, the energies in eV corresponding to this xanes data
%
%ElementName, a string of the element of the absorbing atom of interest
%       e.g. 'Co' or 'O'
%
%E_ROI_bounds, a 2-element vector specifying the E_ROI high and low values
%       ("Energy region of interest" for transforming and plotting)
%       default: 100eV above and below input E_given if not given or <2 elements
%
%ifplot, a boolean: 1 if plotting, 0 if not
%       default: does not plot if input is empty
%
%subtractZ, a boolean: 1 if you want to subtract Z, 0 if you do not
%       default: does NOT subtract Z (the atomic number) from ref input
%
%
%OUTPUTS
%ref, the f'  values for the input energy range E_given
%
%if additional outputs desired, currently gives, in order:
%       ref_full, f' over all E, 
%               transformed over E_ROI and spliced from tables elsewhere 
%       all E (the energy values corresponding to these, see below)
%       imf_full, kk transform of ref_full (f' over all E)
%
%The steps in energy are decided from the input data. Before the input data
%are transformed they are interpolated to an evenly spaced energy grid. The
%spacing on this grid is either 0.5eV steps or the smallest energy
%step in the input data, whichever is larger.
%
%
%DEPENDENCIES
%kkrebook2z, a modification of Valerio Lucarini's software 
%       operates over a specified energy region of interest
%       vectorized for speed and with input range functionality
%elements, an open source matlab function that retrieves Z from ElementName
%       mathworks.com/matlabcentral/fileexchange/9921-chemical-elements
%NIST tables, add new fp and fpp .txt files for each element of interest
%       current elements include Al,Co,Fe,La,Mn,O,Ta...
%
%
%
%This software is distributed under the GNU licence agreement
%by Zachary Porter
%last update: 25 Aug 2015
%email: zdp7@cornell.edu
%SLAC National Accelerator Lab
%SSRL MSD Hard X-Rays Department
%Menlo Park, California, U.S.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Inputs

minstep = 0.5;  %0.5eV is smallest E step size, you may want to change this. 
                %It is the only completely arbitrary choice in the program.

%determine the min. step size (greatest common divisor of the steps)
steps = diff(sort(E_given)); steps = steps(steps>0);
if(length(steps)<4)
    interpmethod = 'linear'; else interpmethod = 'spline';
end
if ~isempty(steps)
    steps = round(1e6*steps); %for the gcd function (divide by this later)
    E_step = min(gcd(steps, min(steps)));
    E_step = max(E_step/(1e6), minstep); %enforce a minimum step size of .5eV
else  %if only one E value
    E_step = 5; %iterate in steps of 5eV
end

%adjust inputs given new step size
ref_in = ref;
E_corrected = E_given(1) : E_step : E_given(end);
ref = interp1(E_given,ref,E_corrected,interpmethod);

%find the energy range
E_min  = min(E_corrected); E_max = max(E_corrected);

E_low  = fliplr(E_min:-E_step:20);  % 20 eV is min sensible value in tables
E_low  = E_low(1:end-1);%remove overlapping value

E_high = E_max:E_step:4e5;          %400keV is max value included in tables
E_high = E_high(2:end); %remove overlapping value

E = [E_low, E_corrected, E_high];

%Energy region of interest
if(exist('E_ROI_bounds', 'var') && length(E_ROI_bounds)>1)
    if(E_ROI_bounds(1)>E_given(1) || E_ROI_bounds(end)<E_given(end))
        error('Select bounds beyond the input energy range');
    end
    E_ROI = E_ROI_bounds(1):E_step:E_ROI_bounds(end);
else
    E_ROI = (E_min-100): E_step :(E_max+100);
end

%adjust element name to 1st letter capital, 2nd lowercase
switch length(ElementName)
    case 1
        ElementName = upper(ElementName);
    case {2, 3}
        ElementName = [upper(ElementName(1)), lower(ElementName(2:end))];
    otherwise
        error('ElementName is IUPAC chemical symbol (1-2 characters)');
end


%% load element data from NIST
try
    fp = load([ElementName,' fp.txt']).';
    fpp = load([ElementName,' fpp.txt']).'; %given is negative
catch
    error(['Need to add element fp and fpp files in folder on path',...
            'and put in same format as existing files']);
end

Z = elements('Symbol', ElementName, 'atomic_number'); Z = double(Z);

fp_low = interp1(fp(1,:).*1e3, fp(2,:), E_low) -Z;
fp_high = interp1(fp(1,:).*1e3, fp(2,:), E_high) -Z;

%subtract Z from real[f] to get f', if not done already
ref = ref-Z;
if (nargin<6 || ~subtractZ) 
    ref = ref+Z;
end

ref_full = [fp_low, ref, fp_high];

%% KK transform into f'
imf_ROI = kkimbook2z(E, ref_full, 0, E_ROI(1), E_ROI(end));
imf = interp1(E, imf_ROI, E_given, interpmethod);
imf_ROI = interp1(E, imf_ROI, E_ROI, interpmethod);


%% Plotting
if (nargin>4 && ifplot)
    
    fpp = interp1(fpp(1,:)*1e3, -1.*fpp(2,:), E_ROI);%given is negative
    figure;
    plot(E_ROI, fpp , 'k'); hold on;
    plot(E_ROI, imf_ROI, 'g', 'LineWidth',2);
    title(['Imag. K-K Calculation for f", ', ElementName, ' Edge']);
    xlabel('Energy (eV)'); ylabel('f"');
    legend({'table', 'KK of f'''}, 'Location', 'southeast');
    
    fp = interp1(fp(1,:)*1e3, fp(2,:), E_ROI)-Z;
    figure;
    plot(E_ROI, fp, 'k'); hold on;
    plot(E_given, ref_in, 'g', 'LineWidth',2);
    title(['Input f'' data, ', ElementName, ' Edge']);
    xlabel('Energy (eV)'); ylabel('f''');
    legend({'table', 'input'}, 'Location', 'southeast');
end

%% Variable output handling
nout = max(nargout,1) - 1;
varargout = cell(1,nout);
for k =1:nout
    switch k
        case 1
            varargout{k} = ref_full;
        case 2
            varargout{k} = E;
        case 3
            imf_ROI = kkimbook2z(E, ref_full, 0);
            varargout{k} = imf_ROI;
        otherwise
            varargout{k} = [];
    end
end

end
