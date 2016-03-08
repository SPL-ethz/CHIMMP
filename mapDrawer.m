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
% This function calculates the morphology map based on the boundary
% conditions found using morphologyDomainExplorer. It returns the map and
% suggests a colormap via output Cmat.
%
% Last modified:
% - 12. Jan 2016, DO: Initial creation
%
%
% Input arguments:
% - A:                      (m x 3)-matrix containing all unit normal vectors of the m possible facets, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
% - M:                      (m x n)-matrix needed for correct grouping of facets into sets that grow with the same speed, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
%
% Output argument:
% - rvec:                   Range of aspect ratios that was investigated / will be the axes in the map
% - Map:                    The actual morphology domain map, an (n-1)-dimensional array. Each value corresponds to a separate (sub-)morphology
% - Cmat:                   Color map associated to morphology map used for plotting purposes
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

function [rvec,Map,Cmat,outFlag] = mapDrawer(A,M)


%% Preparation / Parameters
% Default output flag is 1
outFlag = 1;

if size(M,2) == 2
    outFlag = 0; % there are only true independently growing facet groups... trick a bit 
    A = [A;A(logical(M(:,2)),:)];
    M = [M zeros(size(M,1),1);repmat([0 0 1],[size(M(logical(M(:,2)),:),1),1])];
    
elseif size(M,2) == 1
    outFlag = -3; % there is no morphology map because there is only one independently growing facet
    rvec = [];
    Map = [];
    Cmat = [];
    return;
end

% Prepare A to reduce numerical issues
A(abs(A)<1e-5) = 0;
A = normc(A')';

[~,n] = size(M); % extract info regarding number of facets
[rows,Igroup] = find(M); % store which facet belongs to which group
[~,Isort] = sort(rows);
Igroup = Igroup(Isort);



% investigated aspect ratios for morphology domain
sizeArr = zeros(1,n-1);
rvec = cell(1,n-1);

for i = 1:n-1
        rvec{i} = logspace(log10(exp(-4)),log10(exp(4)),(n<3)*500+(n>=3)*101);
        sizeArr(i) = (n<3)*500+(n>=3)*101;
end
% number of grid points used to draw the map


% quickly check how many morphologies will be present.
% This is a simple sanity check: nMorphsTest must be equal to the number of
% morphology domains found by the code (nMorphs); this condition is necessary, but
% not sufficient for the morphology map to be correct
nMorphsTest = 0;
for i = 1:n
    famCombs = nchoosek(1:n,i);
    for j = 1:size(famCombs,1)
        Asub = A(ismember(Igroup,famCombs(j,:)),:);
        Psub = Polyhedron(Asub,ones(size(Asub,1),1));
        if Psub.isBounded
            nMorphsTest = nMorphsTest + 1;
        end
    end
end

if ~nMorphsTest
   warning('ADDICT:mapDrawer:noBoundedShape',...
       'The morphology provided cannot form a closed polyhedron. The reason is most likely that some fast growing facet was discarded for the calculation.');
   outFlag = -2;
   return;
end

% Create nD arrays of rvec
R       =   cell(1,n-1);
[R{:}]  =   ndgrid(rvec{:});

for i = 1:length(R)
    % ndgrid has a different behavior than meshgrid, we permute so that
    % changes in horizontal axis (column of matrix) correspond to changes
    % in x1
    R{i} = permute(R{i},[2 1 3]);
end

% Compute boundary matrix
coneData    =   morphologyDomainExplorer(A,M);
Bloc        =   coneData.B; % if the edge classification would work perfectly this could be replaced with conData.BfaceGone... alas, it sometimes fails, hence we go through all conditions, even those that should be of newEdge type

% Check which points in the map satisfy which conditions in B
conditionCheck = [];
for i = 1:size(Bloc,1)
    
    mLoc = Bloc(i,1);
    
    for j = 2:n
        mLoc    =   mLoc + Bloc(i,j)*R{j-1}*(j<=3)+Bloc(i,j).*R{j-1}.*(j>3).*R{2};
    end
    
    mLocI           =   false(size(R{1}));
    mLocI(mLoc<=0)  =   true;
    conditionCheck     =   cat(n,conditionCheck,mLocI);
    
end

% reshape array into matrix where each row corresponds to one coordinate
% and columns correspond to whether or not associated condition is
% fulfilled
coordCondInfo    = reshape(conditionCheck,[numel(R{1}) size(Bloc,1)]);

% write down to which of the x unique combinations of condition-fulfillment
% each coordinate belongs
[~,~,coordCondClass]    = unique(coordCondInfo,'rows');

% reshape back into array format
protoMorphMap    =   reshape(coordCondClass,sizeArr);

% Tell us for each candidate morphology domain which facets are present
[realFacets]    =   determineRealFacets(protoMorphMap,R,A,M);

% Tell us for each candidate morphology to which unique submorphology it belongs
% Submorphology: (different sets of facets presents; not necessarily different _groups_ of facets (as grouped by M))
[uniqueRealFacets,~,realFacetClass] = unique(realFacets','rows');
nSubMorph = size(uniqueRealFacets,1); % number of submorphologies

% Now we can create the map of all submorphologies 
subMorphMap = zeros(size(protoMorphMap));
for i = 1:length(realFacetClass)
    subMorphMap(protoMorphMap==i) = realFacetClass(i);
end

% Tell us which groups are present for each submorphology
dummy = (uniqueRealFacets'.*repmat(Igroup,[1 nSubMorph]))';
labelGroups = false(nSubMorph,n);
for i = 1:nSubMorph
    labelGroups(i,unique(dummy(i,dummy(i,:)~=0))) = true;
end

% Determine final class for each submorphology (each class = separate morphology domain)
[~,~,morphClass] = unique(labelGroups,'rows');
nMorphs = max(morphClass); % number of classes

% perform sanity test
if nMorphsTest ~= nMorphs
    warning('mapDrawer:sanityFail',...
        'The built-in sanity check of the morphology map calculator failed. Number of morphology domains was not equal to the expected value.');
    outFlag = -1;
end
% Calculate actual map
Cmat = zeros([size(subMorphMap) 3]); % color matrix
clrs = [111 183 214;
    249 140 182;
    133 202 93;
    255 255 176;
    255 255 255;
    191 213 232;
    255 237 81;
    192 186 153;
    252 169 133; 
    163 137 193;
    145 210 144;
    165 137 193;
    221 212 213;
    179 226 221;
    154 206 223;
    191 213 232]/255;

Map = zeros(size(subMorphMap)); % memory pre-allocation
newSize = ones(1,n);newSize(n) = 3;
for i = 1:nMorphs
    I = find(morphClass == i); % find submorphologies belonging to this class

       helpC = ismember(subMorphMap,I);

       Cmatloc                  = repmat(helpC,newSize).*repmat(reshape(clrs(i,:),newSize),[size(subMorphMap) 1]); % create local color matrix
       Cmat                     = Cmat+Cmatloc; % add to overall color matrix
       Map(ismember(subMorphMap,I))   = i; 
   
end

