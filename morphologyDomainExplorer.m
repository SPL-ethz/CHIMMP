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
% This code explores the morphology domain and returns joint boundary 
% matrix and the joint morphology cone using an algorithm that
% was proposed by Borchert et al in Chem. Eng. Sci. 84, 2012, 85-99. All
% equations that are referenced in the comments refer to that paper. 
% The algorithm essentially solves for conditions where edges/facets
% disappear and explores the morphology domain in a smart way.
% Nevertheless, the  code uses several heuristics and is not guaranteed to 
% return the complete morphology domain, though it has been successfully 
% tested for a variety of cases.
%
% Last modified:
% - 11. Jan 2016, DO: Initial creation
%
%
% Input arguments:
% - A:                      (m x 3)-matrix containing all unit normal vectors of the m possible facets, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
% - M:                      (m x n)-matrix needed for correct grouping of facets into sets that grow with the same speed, needed for Polyhedron description P = {x \in R3 | Ax <= ML})
%
% Output argument:
% - coneData:               Structure describing the morphology cone
%   * morphology            Structure containing information about submorphologies 
%   * Cone                  Joint morphology cone
%   * B                     Joint boundary matrix
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

function [coneData]=morphologyDomainExplorer(A,M)

%% Preparation
A(abs(A)<eps) = 0; % remove numerical artifacts sometimes present in A

[m,n] = size(M); % get m (number of potential facets) and n (number of independently growing groups)

Nijk = @(ijk) [A(ijk(1),:);A(ijk(2),:);A(ijk(3),:)]; % subset of normal vectors (eq. 4 Borchert 2012)

Ltest = ones(n,1); % matrix containing all morphology vectors that shall be tested, starts with unity-vector (all facets are present)

morphology = struct; % structure containing info regarding found morphologies
nMorphs = 0; % counter

while ~isempty(Ltest)
    % code will run until there are no more morphology vectors to be tested 
    % (i.e., Ltest is empty)
    
    L = Ltest(:,1); % the L vector we currently look at
    
    % Compute set of appearing vertices
    calcVertices;

    % Compute set of appearing edges
    calcEdges;
    
    %check whether polytope is new
    newone = true;
    if nMorphs>0
        for i = 1:length(morphology)
            if isequal(morphology(i).E.BorI,E.BorI) && isequal(morphology(i).V.I,V.I) % does the new polytope have the same list of edges and vertices as a previous one?
               newone = false; % this is not a new morphology
            end
        end
    end
    
    if newone % this is in fact a new morphology
        P(nMorphs+1) = Polyhedron(A,M*L); % create the polyhedron
        
        %% Compute Boundary matrix for cone
        % pre-allocation
        B = zeros(E.N,n);BnewEdge = zeros(E.newEdge.N,n);BfaceGone = zeros(E.faceGone.N,n);
        
        % reduce conditions to symmetrical case (eq. 32)
        for i=1:E.N
            B(i,:) = (M'*bfun(E.BorI(i,:)))';
        end

        for i = 1:E.newEdge.N  
            BnewEdge(i,:) = (M'*bfun(E.newEdge.BorI(i,:)))';
        end

        for i = 1:E.faceGone.N
            BfaceGone(i,:) = (M'*bfun(E.faceGone.BorI(i,:)))'; 
        end
        
        % normalization and removing numerical artifacts
        BfaceGone = normc(BfaceGone')';BfaceGone(abs(BfaceGone)<eps) = 0;BfaceGone = unique(BfaceGone,'rows');
        B = normc(B')';B(abs(B)<eps) = 0;B = unique(B,'rows');
        BnewEdge = normc(BnewEdge')';BnewEdge(abs(BnewEdge)<eps) = 0;BnewEdge = unique(BnewEdge,'rows');
        
        % save the newly found submorphology
        morphology(nMorphs+1).B = B;
        morphology(nMorphs+1).BnewEdge = BnewEdge;
        morphology(nMorphs+1).BfaceGone = BfaceGone;
        morphology(nMorphs+1).E = E;
        morphology(nMorphs+1).V = V;
        morphology(nMorphs+1).P = P(nMorphs+1);
        
        
        %% find new guesses for h (eq.26)
        if n > 1 % only makes sense if number of independently growing groups is larger than 1
            for j = 1:size(BnewEdge,1) % for all newEdge boundaries
                
                % find a h on a boundary
                % the problem is that for most cases there is no unique
                % choice... this algo will try to find ONE possible point
                % and it uses a heuristic to deal with special cases
                
                % make a guess and assume there is a point that lies on the boundary whose
                % coordinates are all unity except for one element
                xlast = -1;itried = 1;
                hnew = ones(n,1);                
                while (abs(BnewEdge(j,:)*hnew)>1e-5 || isnan(xlast)) && itried <=n  
                    xlast = fminbnd(@(x) abs(BnewEdge(j,:)*circshift([x; ones(n-1,1);],[itried-1,0])),0,1e5); 
                    hnew = ones(n,1);
                    hnew(itried) = xlast;
                    itried = itried +1;
                end

                    
                if any(isnan(hnew)) ||  abs(BnewEdge(j,:)*hnew)>1e-5
                    % we did not manage to find a suitable point, discard
                    hnew = [];    
                    % DRO: this has so far not created any issues, but could this be a source of problems?
                    
                else % we have found a valid point on the border
                    hnew = hnew/norm(hnew); % normalize
                    hnew = hnew(:) + 0.05*BnewEdge(j,:)'; % move a bit outside the border
                    hnew(hnew<0) = 0;
                    hnew0 = hnew;
                
                    
                    % If an element of BnewEdge is zero, it means that this
                    % facet plays no role in the transition to another
                    % morphology. However, one should check out the
                    % morphologies obtained by varying that distance
                    % anyway. We have missed regions in the past due to
                    % unfortunate choices of the facet distance.
                    if any(BnewEdge(j,:)==0) 
                        for ii = 1:length(BnewEdge(j,:))
                            if BnewEdge(j,ii) == 0
                               haddnew = repmat(hnew0(:),[1 5]);
                               haddnew(ii,:) = hnew(ii)*exp(linspace(-3,3,5));
                               hnew = [hnew haddnew];
 
                            end
                        end
                    end
                    
                end
                
                Ltest = [Ltest hnew]; % add to our list of to be investigated hs
                
            end % for all newEdge boundaries
        end % if n > 1
        
        nMorphs = nMorphs+1; % increase the morphology counter by one
        
    end % if we have found a new morphology (newone = true)
    
    Ltest(:,1) = []; % we can now delete the tested h from the list
    Ltest = unique(Ltest','rows')'; % remove multiple copies of to-be-tested vectors
    
end
    
%% Joint Morphology cone
B = [];
for i = 1:length(morphology)
   B = [B;morphology(i).B];
end
B(abs(B)<eps) = 0;B = unique(B,'rows');

% save joint cone info
coneData.B = B;
coneData.morphology = morphology;

    function calcVertices
        
        V = struct('loc',zeros(0,3),'I',[],'N',0); % reset structure V
        % loc: contains 3D coordinates of all N vertices (Nx3 matrix)
        % I: indices of facets that form a given vertex (Nx3 matrix)
        % N: total number of vertices
        
        % Generate list of all facet combinations that include 3 facets (at
        % least three facets are needed to create a vertex)
        facetCombs  = nchoosek(1:m,3); 
        
        % normal distances for all m facets
        h           = M*L;
                
        for i = 1:size(facetCombs,1)

           Nijk_loc     = Nijk(facetCombs(i,:)); % find subset of normal vectors for this combination of facets
           
            if abs(det(Nijk_loc))>10*eps 
                % check determinant to see if orientations are independent
                
                % calculate locations of vertices (eq. 6)
                hijk     = [h(facetCombs(i,1),:);h(facetCombs(i,2),:);h(facetCombs(i,3),:)]; 
                vijk     = Nijk_loc\hijk; 

                % add newly found vertex to list if it is not virtual
                if all(A*vijk<=M*L+1e-5)
                    V.loc   = [V.loc; vijk'];
                    V.I     = [V.I; facetCombs(i,:)];
                end
                
            end % if determinant not zero
        end % for all combinations of facets
        
        V.loc(abs(V.loc)<eps) = 0; % get rid of numerical issues
        V.N = length(V.loc(:,1)); 

    end

    function calcEdges
        
        Btemplate = struct('N',0,'BorI',zeros(0,4),'VI',zeros(0,2));
        E = struct('BorI',zeros(0,4),'N',0,'VI',zeros(0,2),'faceGone',Btemplate,'newEdge',Btemplate); % reset structure E
        % BorI: indices of facets that form ith edge
        % N: total number of edges
        % VI: indices of vertices that form ith edge
        % newEdge: boundary vectors for boundaries where new edges appear
        % faceGone: boundary vectors for boundaries where facets disappear
        
        % find all vertices that share a facet -> edge
        vertCombs = nchoosek(1:V.N,2); % determine all possible combinations of vertices
        
        for i = 1:size(vertCombs,1)
            commonFacets = myIntersect(V.I(vertCombs(i,1),:),V.I(vertCombs(i,2),:)); % find indices of facets that a pair of vertices has in common
            
            if length(commonFacets) > 1 && norm(V.loc(vertCombs(i,1),:)-V.loc(vertCombs(i,2),:))>eps% if they have at least two facets in common there is an edge between the two vertices
  
               E.N = E.N+1; % congratulations, you have found a new edge
               differentFacets = mySetXOr(V.I(vertCombs(i,1),:),V.I(vertCombs(i,2),:),commonFacets); % which facets are not part of the edge. k and l indices in borchert nomenklatur

               % make sure Nijk is a righthanded system and Nijl a left-handed system 
               if det(Nijk([commonFacets differentFacets(1)]))<0 % determinant negative -> left-handed system
                   differentFacets = circshift(differentFacets(:)',[0 1]);
               end

               if isempty(find(V.I(vertCombs(i,1),:)==differentFacets(1), 1)) % ensure right order
                  vertCombs(i,:) = circshift(vertCombs(i,:)',[0 1]) ;
               end

               E.BorI   = [E.BorI; commonFacets differentFacets]; % for each edge: write down the involved facets
               E.VI     = [E.VI;vertCombs(i,1),vertCombs(i,2);];
               
               % classify the latest edge (eq. 10) / checks for determinant
               % because otherwise matrix is too close to being singular
               % note that this test understimates the number of
               % transitions that lead to facet disapperance (the positive
               % linear span condition seems to be sufficient but not
               % necessary for the symmetrical case?!). Regardless, this
               % just slows down the exploration a bit, but not
               % dramatically.

               if (abs(det((Nijk(E.BorI(end,2:4))')))<10*eps || all((Nijk(E.BorI(end,2:4))')\(A(E.BorI(end,1),:)')>=0)) || (abs(det((Nijk(E.BorI(end,[1 3 4]))')))<10*eps || all((Nijk(E.BorI(end,[1 3 4]))')\(A(E.BorI(end,2),:)')>=0))
                  E.faceGone.N      = E.faceGone.N + 1;
                  E.faceGone.VI     = E.VI(end,:);
                  E.faceGone.BorI   = [E.faceGone.BorI;E.BorI(end,:)];
               else
                  E.newEdge.N       = E.newEdge.N + 1;
                  E.newEdge.VI      = E.VI(end,:);
                  E.newEdge.BorI    = [E.newEdge.BorI; E.BorI(end,:)];
               end

                   
           end %distance>eps
            
        end % for all pairs of vertices
        
    end

    function [bijkl] = bfun(ijkl)
        
        bijkl = zeros(m,1);
        
        bijkl(ijkl(1)) = -det(Nijk(ijkl([2 3 4])));
        bijkl(ijkl(2)) =  det(Nijk(ijkl([1 3 4])));
        bijkl(ijkl(3)) = -det(Nijk(ijkl([1 2 4])));
        bijkl(ijkl(4)) =  det(Nijk(ijkl([1 2 3])));

        bijkl = -bijkl/norm(bijkl);

    end

end

function C = myIntersect(A,B)
% this function should be faster in calculating the intersection than
% intersect in the case we are interested in

if ~isempty(A)&&~isempty(B)
    P       = false(1, max(max(A),max(B)) ) ;
    P(A)    = true;
    C   = B(P(B));
else
    C   = [];
end

end

function C = mySetXOr(A,B,Cintersect)

C = unique([A B]);
C = C(~ismember(C,Cintersect));

end