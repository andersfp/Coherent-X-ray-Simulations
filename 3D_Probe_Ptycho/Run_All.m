% Initialization
clear;
close all;
clc;


%% Run script
% Set the path names
p = {...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\Lens_Aberrations_Dislocations_Phase\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\Lens_Aberrations_No_Phase\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\No_Aberrations_Dislocations_Phase\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\No_Aberrations_No_Phase\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\NA_2x\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\NA_4x\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\NA_8x\',...
    %'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\Partial_Coherence\',...
    'C:\Users\anfils\OneDrive\DTU\PostDoc\Projects\3D_Probe_Ptycho\Partial_Coherence_2\',...
    };
np = length(p);

% Set the file names
f = {'Generate_Phantom_2.m','Generate_Exit_Fields.m','Forward_Propagation.m','Reconstruction.m'};
nf = length(f);

% Combine path and file names
lst = cell(nf,np);
for i = 1:nf
    for j = 1:np
        lst(i,j) = {[p{j} f{i}]};
    end
end
lst = lst(:);

% Add the projection and probe only script
%lst{end+1} = [pwd '\Projection_Probe_Only.m'];

% Make command list
cmd = cellfun(@(x) ['run(''' x ''');'],lst,'UniformOutput',false);

% Combine all commands
cmd2 = '';
for i = 1:(nf*np)
    cmd2 = [cmd2 cmd{i}];
end

% Execute command list
eval(cmd2);

