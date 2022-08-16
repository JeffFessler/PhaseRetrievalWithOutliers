classdef (ConstructOnLoad=true) myDFT
    %myDFT This is a stripped-down version of Gdft from the IRT that does
    %what I need for simulations, without some of the more complicated
    %fatrix support. Stripping out a lot of the fatrix functionality
    %actually greatly speeds up computation when we're computing lots of
    %DFT's for small signals (like in the Monte Carlo simulations). 
    %
    %Gdft is part of the Image Reconstruction Toolbox, copyright Jeffrey
    %Fessler and the University of Michigan. The toolbox may be downloaded
    %at http://web.eecs.umich.edu/~fessler/code/index.html. 
    
    properties (SetAccess = protected)
        N = 0
        FFTsize = 0
        samp = []
        adjoint = false
        scale = 1
    end
    
    methods
        function obj = myDFT(N,FFTsize,samp)
            
            if exist('N','var') && ~isempty(N), obj.N = N(1); obj.FFTsize = obj.N; end
            if exist('FFTsize','var') && ~isempty(FFTsize), obj.FFTsize = FFTsize(1); end
            if exist('samp','var') && ~isempty(samp), obj.samp = samp(:); end
            
            if ~isempty(obj.samp) && length(obj.samp) ~= obj.FFTsize, error('# of samples must match FFTsize!'); end
            
        end
        
        function y = ctranspose(obj)

            y = obj;
            y.adjoint = ~y.adjoint;
        end
        
        function varargout = size(obj,dim)
            if isempty(obj.samp)
                sz = [obj.FFTsize,obj.N];
            else
                sz = [sum(obj.samp),obj.N];
            end
            if ~exist('dim','var')
                nout = nargout;
                sz(end+1:nout) = 1;
                if nout > 1
                    varargout(1:nout-1) = num2cell(sz(1:nout-1));
                end
                varargout{nout} = sz(nout:end);
            else
                if dim <= length(sz)
                    varargout{1} = sz(dim);
                else
                    varargout{1} = 1;
                end
            end
        end
        
        function y = mtimes(A,B)
            
            if isa(A,'myDFT')
                obj = A; x = B; pre = false;
            else
                obj = B; x = A; pre = true;
            end
            
            if ~isnumeric(x), error('Object*Object not supported!'); end
            
            if isscalar(x) % 
                y = obj;
                y.scale = y.scale * x;
                return;
            end
            
            if pre, error('x*Object only supported for scalar'); end

            sz_x = size(x);
            
            if obj.adjoint % y = DFT'*x
                if isempty(obj.samp)
                    sz_L = sz_x(2:end); L = prod(sz_L);
                    y = reshape(x,[obj.FFTsize,L]);
                else
                    sz_L = sz_x(2:end); L = prod(sz_L);
                    y = zeros([obj.FFTsize,L],class(x));
                    y(obj.samp,:) = reshape(x,[],L);
                end
                y = ifft(y,obj.FFTsize,1);
                y = obj.FFTsize.*reshape(y(1:obj.N,:),[obj.N,sz_L,1]);
                y = conj(obj.scale).*y;
            else % y = DFT*x
                sz_L = sz_x(2:end); L = prod(sz_L);
                y = fft(x,obj.FFTsize,1);
                if isempty(obj.samp)
                    y = reshape(y,[obj.FFTsize,sz_L]);
                else
                    y = reshape(y,[obj.FFTsize,L]);
                    y = y(obj.samp,:);
                    y = reshape(y,[size(y,1),sz_L]);
                end
                y = obj.scale.*y;
            end
            
        end
        
    end
    
end

