function [pFace,MilCol,Pmin,clrs] = colorCrystal(Pbase,P,M,Mil,clrs)

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

Pmin = Polyhedron(P.A,P.b);
Pmin.minHRep;

nGroups = size(M,2);

if nargin < 5
    clrs    = distinguishable_colors(nGroups,{'k' 'w'});
end
MilCol  = zeros(nGroups,6);

if nargout>1
    for j = 1:size(Pbase.A,1)
        [Icolor]= find(M(j,:));

        if all(MilCol(Icolor,:)==0) || sum(MilCol(Icolor,1:3)<0)>sum(Mil(j,1:3)<0) || (sum(MilCol(Icolor,1:3)<0)==sum(Mil(j,1:3)<0) && sum(MilCol(Icolor,1:3).*[3 2 1])>sum(Mil(j,1:3).*[3 2 1]))
           MilCol(Icolor,:) = [Mil(j,:) clrs(Icolor,:)];
        end
    end
end

pFace = [];
for j = 1:length(Pmin.A(:,1))
    [~,I] = ismember(Pmin.A(j,:),Pbase.A,'rows');
    [~,Icolor]= find(M(I,:));

    pFace = [pFace; clrs(Icolor,:)];
end