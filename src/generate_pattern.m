function [image,dictionary,coeffs,specs] = generate_pattern(file,pattern,varargin)
% [image,dictionary,coeffs,specs] = generate_pattern(file,pattern,key1,value1,key2,value2,...)
%
% Inputs:
%  file - file containing points for pattern centers in image
%  pattern - pattern to place in image and dictionary (like a disc)
%  keys, values - see "opts" structure below
%
% Outputs:
%  image - true image
%  dictionary - set of atoms corresponding to shifts of pattern
%  coeffs - true sparse coefficients (image = dictionary*coeffs)
%  specs - additional info about pattern and image
%
% This function constructs an image corresponding to arranging the
% specified pattern according to the center points listed in the file.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

opts.border = 1; % empty border around image; use [x y] for different borders
opts.sz_pow2 = true; % make image size a power of 2
opts.sz_square = false; % make image size square

opts = vararg_pair(opts,varargin);

%% read pattern info from file
fid = fopen(file,'r');
centers = fscanf(fid,'%g',[2,Inf]).';
fclose(fid);

%% scale by pattern size
sz_pattern = size(pattern);

centers = round(bsxfun(@times,centers,sz_pattern));

%% calculate image size
sz_image = ceil(2*max(abs(centers),[],1)+sz_pattern+2.*opts.border);
if opts.sz_pow2
    sz_image = 2.^nextpow2(sz_image);
end
if opts.sz_square
    sz_image = max(sz_image).*ones(1,2);
end
coeffs = zeros(sz_image);

centers = bsxfun(@plus,centers,floor(sz_image/2))+1;

%% put in coeffs
centers = sub2ind(sz_image,centers(:,1),centers(:,2));
coeffs(centers) = 1;

%% construct dictionary (uses convolution) and image
pattern_fft = zeros(sz_image);
pattern_fft(ceil((sz_image(1)-sz_pattern(1))/2) + (1:sz_pattern(1)),ceil((sz_image(2)-sz_pattern(2))/2) + (1:sz_pattern(2))) = pattern;
pattern_fft = fftshift(fft2(ifftshift(pattern_fft)));
pattern_scale = max(abs(pattern_fft(:)));

dictionary = fatrix2('odim',sz_image,'idim',sz_image,'arg',struct('sz_image',sz_image,'pattern_fft',pattern_fft./pattern_scale),'forw',@generate_pattern_dictionary_forw,'back',@generate_pattern_dictionary_back);
coeffs = pattern_scale.*coeffs;

image = abs(dictionary*coeffs);

%% save specs
specs = struct('pattern_size',sz_pattern,'image_size',sz_image,'border',opts.border,'pattern',pattern./pattern_scale,'pattern_fft',pattern_fft./pattern_scale);

end

function y = generate_pattern_dictionary_forw(arg,x)

x = reshape(x,arg.sz_image); % make image
y = fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(x))).*arg.pattern_fft)));

end

function z = generate_pattern_dictionary_back(arg,y)

y = reshape(y,arg.sz_image);
z = fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(y))).*conj(arg.pattern_fft))));

end

