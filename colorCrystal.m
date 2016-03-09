function [pFace,MilCol,Pmin] = colorCrystal(Pbase,P,M,Mil)

Pmin = Polyhedron(P.A,P.b);
Pmin.minHRep;

nGroups = size(M,2);

clrs = distinguishable_colors(nGroups,{'k' 'w'});
MilCol = zeros(nGroups,6);

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