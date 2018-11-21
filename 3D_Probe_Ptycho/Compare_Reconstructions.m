% Initialization
clear;
close all;
clc;

% Save videos?
sav = 0;


%% Load data
% Set the paths
p = {...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_Dislocations_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Lens_Aberrations_No_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_Dislocations_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\No_Aberrations_No_Phase\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_2x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_4x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\NA_8x\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence\',...
    'C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Partial_Coherence_2\',...
    };

% Select data sets
p = p([1 6 8 9]);

% Load the true object
load([p{1} 'Phantom.mat']);
nx = length(x);
ny = length(y);
nz = length(z);

% Load the reconstructions
n = length(p);
R(n) = struct('rho',[]);
for i = 1:n
    R(i) = load([p{i} 'Reconstruction.mat'],'rho');
end

% Load the probe/projection only
P = load('C:\Users\anfils\Documents\Simulation_Results\3D_Probe_Ptycho\Projection_Probe_Only.mat','rho');
P = P.rho;


%% Process the data
% Adjust the phase
obj = obj.*exp(-1i.*angle(obj(ny/2+1,nx/2+1,nz/2+1)));
for i = 1:n
    R(i).rho = R(i).rho.*exp(-1i.*angle(R(i).rho(ny/2+1,nx/2+1,nz/2+1)));
end

% Normalize the reconstructions
P = P./P(ny/2+1,nx/2+1,nz/2+1);
for i = 1:n
    R(i).rho = R(i).rho./(mean(abs(R(i).rho(abs(R(i).rho) > mean(abs(R(i).rho(:)))))) + std(abs(R(i).rho(abs(R(i).rho) > mean(abs(R(i).rho(:)))))));
end


%% Plot the data in 3D
% Axis scale factor
sc = 1e6;

% 3D plot the true object
Slicer(abs(obj),'displayRange',[0 1]);
Slicer(angle(obj),'displayRange',[-pi pi]);

% 3D plot the projection only
Slicer(P,'displayRange',[0 1]);

% 3D plot the reconstructions
for i = 1:n
    Slicer(abs(R(i).rho),'displayRange',[0 1]);
    Slicer(angle(R(i).rho),'displayRange',[-pi pi]);
end

% Make a 3D rendering of the object
figure;
pp = patch(isosurface(sc.*rx,sc.*ry,sc.*rz,abs(obj),0.5));
isonormals(sc.*rx,sc.*ry,sc.*rz,abs(obj),pp);
pp.FaceColor = 'red';
pp.EdgeColor = 'none';
daspect([1 1 1]);
view(3);
camlight(45,45);
camlight(-135,-45);
lighting gouraud;
xlabel('x [\mum]');
ylabel('y [\mum]');
zlabel('z [\mum]');
if sav == 1
    dt = 1;
    rot = [cosd(dt) -sind(dt) 0;sind(dt) cosd(dt) 0;0 0 1].';
    ax = gca;
    ax.CameraViewAngleMode = 'manual';
    ax.XLim = [-2.5 2.5];
    ax.YLim = [-2.5 2.5];
    ax.ZLim = [-2.5 2.5];
    set(gcf,'Color','white','Position',[460 100 1000 1000]);
    nr = round(360./dt);
    mov(nr) = struct('cdata',[],'colormap',[]);
    for i = 1:nr
        ax.CameraPosition = ax.CameraPosition*rot;
        drawnow;
        mov(i) = getframe(gcf);
    end
    vid = VideoWriter('True_Object.mp4','MPEG-4');
    vid.FrameRate = 60;
    open(vid);
    writeVideo(vid,mov);
    close(vid);
end

% Make 3D renderings of the reconstructions
pos = [0 0.5;0.5 0.5;0 0;0.5 0];
figure;
for i = 1:n
    subplot(2,2,i)
    pp = patch(isosurface(sc.*rx,sc.*ry,sc.*rz,abs(R(i).rho),0.5));
    isonormals(sc.*rx,sc.*ry,sc.*rz,abs(R(i).rho),pp);
    pp.FaceColor = 'red';
    pp.EdgeColor = 'none';
    daspect([1 1 1]);
    view(3);
    camlight(45,45);
    camlight(-135,-45);
    lighting gouraud;
    xlabel('x [\mum]');
    ylabel('y [\mum]');
    zlabel('z [\mum]');
    set(gca,'OuterPosition',[pos(i,:) 0.5 0.5]);
end
if sav == 1
    ax = get(gcf,'Children');
    [ax.CameraViewAngleMode] = deal('manual');
    [ax.XLim] = deal([-2.5 2.5]);
    [ax.YLim] = deal([-2.5 2.5]);
    [ax.ZLim] = deal([-2.5 2.5]);
    set(gcf,'Color','white','Position',[460 100 1000 1000]);
    mov(nr) = struct('cdata',[],'colormap',[]);
    for i = 1:nr
        [ax.CameraPosition] = deal(ax(1).CameraPosition*rot);
        drawnow;
        mov(i) = getframe(gcf);
    end
    vid = VideoWriter('Reconstructions_3.mp4','MPEG-4');
    vid.FrameRate = 60;
    open(vid);
    writeVideo(vid,mov);
    close(vid);
end


%% Plot data in 2D
% Set the plotting plane indices
ix = {1:nx,1:nx,129,90};
iz = {129,90,1:nz,1:nz};
xlab = {'x','x','z','z'};
choord = {rx,rx,rz,rz};
tit = {'xy-plane, z = 0 nm','xy-plane, z = -762 nm','zy-plane, x = 0 nm','zy-plane, x = -762 nm'};

% Plot the object amplitude
figure;
for i = 1:4
    subplot(2,2,i);
    imagesc(sc.*squeeze(choord{i}),sc.*ry,squeeze(abs(obj(:,ix{i},iz{i}))),[0 1]);
    axis equal tight;
    xlabel([xlab{i} ' [\mum]']);
    ylabel('y [\mum]');
    title(tit{i});
end
colormap gray;

% Plot the object phase
figure;
for i = 1:4
    subplot(2,2,i);
    imagesc(sc.*squeeze(choord{i}),sc.*ry,squeeze(angle(obj(:,ix{i},iz{i}))),[-pi pi]);
    axis equal tight;
    xlabel([xlab{i} ' [\mum]']);
    ylabel('y [\mum]');
    title(tit{i});
end
colormap hsv;

% Plot the projection amplitude
figure;
for i = 1:4
    subplot(2,2,i);
    imagesc(sc.*squeeze(choord{i}),sc.*ry,squeeze(abs(P(:,ix{i},iz{i}))),[0 1]);
    axis equal tight;
    xlabel([xlab{i} ' [\mum]']);
    ylabel('y [\mum]');
    title(tit{i});
end
colormap gray;

% Plot the reconstruction amplitudes
for j = 1:n
    figure;
    for i = 1:4
        subplot(2,2,i);
        imagesc(sc.*squeeze(choord{i}),sc.*ry,squeeze(abs(R(j).rho(:,ix{i},iz{i}))),[0 1]);
        axis equal tight;
        xlabel([xlab{i} ' [\mum]']);
        ylabel('y [\mum]');
        title(tit{i});
    end
    colormap gray;
end

% Plot the reconstruction phases
for j = 1:n
    figure;
    for i = 1:4
        subplot(2,2,i);
        imagesc(sc.*squeeze(choord{i}),sc.*ry,squeeze(angle(R(j).rho(:,ix{i},iz{i}))),[-pi pi]);
        axis equal tight;
        xlabel([xlab{i} ' [\mum]']);
        ylabel('y [\mum]');
        title(tit{i});
    end
    colormap hsv;
end


%% Other plots
% Plot the strain
figure;
semilogy(sc.*rx(s > 0),s(s > 0));
hold on;
semilogy(sc.*rx(s < 0),-s(s < 0));
xlabel('x [\mum]');
ylabel('Strain');


