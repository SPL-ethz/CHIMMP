%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% UC Santa Barbara, USA
%
% Project:  Novartis (UCSB+ETHZ)
% Year:     2016
% MATLAB:   R2014b, Windows 10
% Authors:  Dave Ochsenbein (DRO)
%
% Purpose:
% This function translates the info received from ADDICT into two matrices
% A and M which are needed for the morphology map calculation. A is
% the set of unit normal vectors of all m possible facets in
% carthesian coordinates; M is the matrix grouping facets of identical
% growth together (allows dimensionality reduction (m->n where n <=m) of scaling vector L).
%
% Last modified:
%
% Input arguments:
% - unitCell:               Structure containing unit cell info
%   * a                     unit cell length a
%   * b                     unit cell length b
%   * c                     unit cell length c
%   * alpha                 unit cell angle alpha [deg]
%   * beta                  unit cell angle beta [deg]
%   * gamma                 unit cell angle gamma [deg]
% - facetInfo:              Array of m structures containing infos regarding individual facets (note that a matrix might be slightly more efficient to read out, but the structure is more versatile for future additions/modifications)
%   * millerIndex           Miller index of facet
%   * G                     Normal growth rate associated to facet (the absolute value here is not of particular importance)
%
% Output arguments:
% - A:                      (m x 3)-matrix containing all unit normal vectors of the m possible facets, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
% - M:                      (m x n)-matrix needed for correct grouping of facets into sets that grow with the same speed, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 David Ochsenbein
% 
% This file is part of CHIMMP.
% 
% CHIMMP is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CHIMMP is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [A,M] = interpretCrystallography(unitCell,facetInfo)


%% Calculate M
m       = length(facetInfo); % number of potential facets

if isfield(facetInfo,'G')
    Gvec = [facetInfo.G]; % write all growth rates into a vector

    [uG,~,In] = unique(Gvec); % find unique values 

    n = length(uG); % number of independently growing groups is equal to unique values in G

    M = zeros(m,n); % pre-allocate memory for mapping matrix
    for i = 1:m
        M(i,In(i)) = 1;
    end

else
    M = [];
end


%% Calculate A
% extract variable names and give them shorter names for readability
alpha   = unitCell.alpha;
beta    = unitCell.beta;
gamma   = unitCell.gamma;
auc     = unitCell.a;
buc     = unitCell.b;
cuc     = unitCell.c;



% Create mapping matrix to map points into cartesian space for (non)cubic
% crystals, assuming c-axis is aligned with z-axis in cartesian space; a-c
% plane lies in x-z plane.
% See the pdf file in this folder explaining the origin of the
% calculations for reference.
delta   = acosd((cosd(gamma)-cosd(beta)*cosd(alpha))/(sind(beta)*sind(alpha)));

Map     = [auc*sind(beta) buc*sind(alpha)*cosd(delta) 0;
            0 buc*sind(alpha)*sind(delta) 0;
            auc*cosd(beta) buc*cosd(alpha) cuc];

% unit cell vectors in cartesian coordinates
a       = Map*[1 0 0]'; 
b       = Map*[0 1 0]';
c       = Map*[0 0 1]';

Vuc     = cross(a,b)'*c; % unit cell volume

% Determine reciprocal lattice vectors
astar   = cross(b,c)/Vuc; 
bstar   = cross(c,a)/Vuc;
cstar   = cross(a,b)/Vuc;

% calculation of unnormalized A -> reciprocal lattice vector in cartesian
% space for all m planes hkl
protoA  = zeros(3,m); % memory pre-allocation
for i = 1:m
    protoA(:,i) = facetInfo(i).millerIndex(1)*astar...
        +facetInfo(i).millerIndex(2)*bstar...
        +facetInfo(i).millerIndex(3)*cstar;
end

A       = normc(protoA)'; % normalization

% Make sure normal of slowest growing facet lies parallel to z-axis
% (assumes first column in M represents slowest facet (ADDICT norm))
R = vrrotvec2mat(vrrotvec([0 0 1],A(find(M(:,1)==1,1,'first'),:)));
A = A*R;


