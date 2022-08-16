function irtpath = setup_IRT(irtpath)
% irtpath = setup_IRT(irtpath)
%
% Inputs:
%  irtpath - path to toolbox (optional)
%
% Outputs:
%  irtpath - path to toolbox, or empty/0 if not found
%
% This function attempts to locate and set up the Image Reconstruction
% Toolbox. If you do not have the toolbox already installed on your MATLAB
% path, this function prompts you to locate the toolbox. Once the toolbox
% is located, it calls setup.m in the toolbox's main directory.
%
% If you do not possess a copy of the Image Reconstruction Toolbox, you may
% download it from http://web.eecs.umich.edu/~fessler/code/index.html.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

im_path = which('im','-all');
[pathstrs,filenames,fileexts] = cellfun(@(p) fileparts(p), im_path,'UniformOutput',false);
goodfiles = strcmp(filenames,'im') & strcmp(fileexts,'.m') & ~cellfun(@isempty,pathstrs);
pathstrs = pathstrs(goodfiles);
pathparts = cellfun(@(p) strsplit(p,filesep),pathstrs,'UniformOutput',false);
goodpaths = cellfun(@(p) strcmp(p{end},'graph'),pathparts);
if any(goodpaths)
    igoodpath = find(goodpaths,1,'first');
    pathparts = pathparts{igoodpath};
    irtpath = fullfile(pathparts{1:end-1});
else
    % need to load IRT
    if ~exist('irtpath','var') || isempty(irtpath), irtpath = ['..' filesep 'irt']; end
    while exist(irtpath,'dir') ~= 7 || exist([irtpath filesep 'setup.m'],'file') ~= 2
        if usejava('awt')
            irtpath = uigetdir('.','Locate Image Reconstruction Toolbox...');
        else
            irtpath = input('Path to IRT: ','s');
        end
        if isequal(irtpath,0)
            return; % later function calls may fail!
        end
    end
    oldp = cd(irtpath);
    setup;
    cd(oldp);
end

end
