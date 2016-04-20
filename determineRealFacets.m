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
% This function determines which of the m potential facets are present in
% all submorphologies of a given morphology map. 
%
% Last modified:
% - 12. Jan 2016, DO: Initial creation
%
%
% Input arguments:
% - morphMap:               Morphology map that should be analyzed.
% - R:                      Cell array of nD grids of aspect ratios associated to morphology map.
% - A:                      (m x 3)-matrix containing all unit normal vectors of the m possible facets, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
% - M:                      (m x n)-matrix needed for correct grouping of facets into sets that grow with the same speed, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
%
% Output argument:
% - facetPresent:           (m,c)-Array showing which facets are present for the c proposed morphologies (unique values in candidate morphology map)
% - Prep:                   Array of representative Polyhedra for all c candidate morphologies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [facetPresent,Prep]=determineRealFacets(morphMap,R,A,M)

[m,n] = size(M); % extract info regarding number of facets
X = cell(1,n-1);
Xavg = cell(1,n-1);


Prep = Polyhedron(1,max(morphMap(:)));
facetPresent = zeros(m,max(morphMap(:)));

for i = 1:max(morphMap(:))
    
    %% look for a representative point in the subdomain i
    
    % find all indices belonging to subdomain
    I = find(morphMap==i);
    [X{:}] = ind2sub(size(R{1}),I);

  
    % use center if possible... looks nicer if you need a
    % 'representative' shape for a given domain
    outsidePoint = false;
    Xavg = round(cellfun(@mean,X));
    if ~ismember(Xavg,[X{1,1:n-1}],'rows');
        outsidePoint = true;
    end
    Xavg = num2cell(Xavg);

    % if the center lies outside the domain, just pick the middle index
    if ~outsidePoint
        Irep = sub2ind(size(R{1}),Xavg{:});
    else
        Irep = I(ceil(length(I)/2));
    end
    
    % create a represntative scaling vector
    if n==3
        Lrep = [1 cellfun(@(c) c(Irep),R)]';
    elseif n ==4
        Lrep = [1 R{1}(Irep) R{2}(Irep) R{2}(Irep)*R{3}(Irep)]';
    end

    % create the corresponding Polyhedron
    Prep(i) = Polyhedron(A,M*Lrep);
    
    % create irredundant A matrix for given polytope
    Amin = Prep(i).minHRep.H; 
    Amin(:,end) = []; % last column not part of A in our terminology

    [facetPresent(:,i)] = ismember(A,Amin,'rows');     
    
end
