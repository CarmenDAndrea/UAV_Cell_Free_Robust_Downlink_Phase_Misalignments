% This function identifies which access points (APs) are linked to each CPU
% in the network. Moreover, it associates the UE k with a master AP and a
% primary CPU.
% Author: Marx Freitas

% INPUT:
% nbrOfCPUs         = the number of CPUs (J) in the newtork
% x_axis_size       = the length in meters of the X-axis of the rectangular coverage area
% y_axis_size       = the length in meters of the Y-axis of the rectangular coverage area
% APpositions_m     = the AP positions in the coverage area in meters
%                     (2D Distance) as a vector of complex numbers
%                     Dim: nbrOfAPs x 1
% nbrOfAPs          = the number of APs (L) in the newtork
% nbrOfUEs          = the number of UEs (K) in the newtork
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% APs_CPUs_id       = the APs linked by fronthaul to each CPU
%                     Dim: nbrOfAPs x nbrOfCPUs
% UEprimaryCPU_ID   = the primary CPU of each UE
%                     Dim: nbrOfUEs x 1

function [APs_CPUs_id, UEprimaryCPU_ID] =...
    MF_Receiving_IDandMasteraAP(nbrOfCPUs, x_axis_size, y_axis_size, APpositions_m,...
    nbrOfAPs, nbrOfUEs, masterAPs)
%% Prepare to store data
APs_CPUs_id = zeros(nbrOfAPs,nbrOfCPUs);
UEprimaryCPU_ID  = zeros(nbrOfUEs,1);
masterOfUE_id = zeros(nbrOfUEs,1);

%=========================================================================%
% ============= ADDING MULTIPLE CPUs IN THE COVERAGE AREA =============== %
%=========================================================================%
% Divide the coverage area into J regions of the same size, where J is the
% number of CPUs, i.e., J = nbrOfCPUs. Thus, j = 1, 2, ..., nbrOfCPUs. Note
% that j is one of the J regions.

% We assume that each region j corresponds to the area controlled by a CPU.
% Hence, we use the same terminology to name the CPUs indexes. Therefore, a
% generic CPU in the network is named CPU j
% -------------------------------------------------------------------------
vertices_CPU = dividingCoverageArea(x_axis_size, y_axis_size, nbrOfCPUs);

%=========================================================================%
% ===================== ASSIGNING APs TO THE CPUs ======================= %
%=========================================================================%
% The access points (APs) inside the region j are controlled by the CPU j
% -------------------------------------------------------------------------
for indexCPU = 1:nbrOfCPUs
    vertices = vertices_CPU{indexCPU}; xv = vertices(:,1); yv = vertices(:,2);
    [APs_CPUs_id(:,indexCPU), ~] = inpolygon(real(APpositions_m),imag(APpositions_m),xv, yv);
end

% Exclude variables to save memory
clear V C bs_ext vertices xv yv

%=========================================================================%
% ============ MASTER AP ASSIGNMENT AND RECEIVING A CPU ID ============== %
%=========================================================================%
% The master AP forwards the identifier (ID) of the primary CPU to the UE
% -------------------------------------------------------------------------
for indexUE = 1:nbrOfUEs

    % Identify the master AP of UE k
    masterOfUE_id(indexUE) = find(masterAPs(:,indexUE) == 1);

    % The master AP assigns the primary CPU identifier (ID) to the UE k
    UEprimaryCPU_ID(indexUE) = find(APs_CPUs_id(masterOfUE_id(indexUE),:) == 1);
end

end

% Local function
function vertices = dividingCoverageArea(x_axis_size, y_axis_size, nbrOfCPUs)
% Vertices of the initial rectangle
verticesIniRec = [0, 0; 0, y_axis_size; x_axis_size, y_axis_size; x_axis_size, 0];

% Computing the number of divisions desired on each axis (x and y)
max1 = floor(sqrt(nbrOfCPUs)); % The first number cannot be greater than the square root of nbrOfCPUs

for nbrOfDivisions_smallest = max1:-1:1
    if rem(nbrOfCPUs, nbrOfDivisions_smallest) == 0
        nbrOfDivisions_largest = nbrOfCPUs / nbrOfDivisions_smallest;
        break;
    end
end

% Atribui o maior número de divisões ao eixo de maior comprimento
if x_axis_size > y_axis_size
    nbrOfDivisions_x = nbrOfDivisions_largest;
    nbrOfDivisions_y = nbrOfDivisions_smallest;
else
    nbrOfDivisions_y = nbrOfDivisions_largest;
    nbrOfDivisions_x = nbrOfDivisions_smallest;
end

% Computing the size of rectangular sub areas
areaSize_x = (verticesIniRec(3, 1) - verticesIniRec(1, 1)) / nbrOfDivisions_x;
areaSize_y = (verticesIniRec(2, 2) - verticesIniRec(1, 2)) / nbrOfDivisions_y;

% Initialize vector to store sub area vertices
vertices = cell(nbrOfCPUs);

% Loop para calcular os vértices das áreas
for i = 1:nbrOfDivisions_x
    for j = 1:nbrOfDivisions_y
        % Calcular vértices da área atual
        x1 = verticesIniRec(1, 1) + (i - 1) * areaSize_x;
        y1 = verticesIniRec(1, 2) + (j - 1) * areaSize_y;
        x2 = x1 + areaSize_x;
        y2 = y1 + areaSize_y;

        % Store the vertices in the matrix
        vertices((i - 1) * nbrOfDivisions_y + j) = {[x1, y1; x1, y2; x2, y2; x2, y1]};
    end
end
end