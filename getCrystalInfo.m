function [crystal,I] = getCrystalInfo(crystalType)
% This function loads the CrystalData table and finds the data stored for
% the crystal type you're looking for.

if ~ischar(crystalType) && ~iscell(crystalType)
        
        error('getCrystalInfo:crystalTypeNotAString',...
            'Variable crystalType must be a string or a cell object (if cell first element is chosen)');
elseif iscell(crystalType)
    
    crystalType = crystalType{1};
    
end
    

try
    load('CrystalData')
catch
    error('getCrystalInfo:failtoload',...
        'Couldn''t find CrystalData.mat. Make sure it is in your search path');
end

I = false(length(Crystals),1);
for i = 1:length(Crystals)
   I(i) = ~isempty(find(strcmpi(crystalType,Crystals(i).heteronyms), 1));
end

if sum(I)>1
    error('getCrystalInfo:multiplyDefinedHeteronyms',...
        'The CrystalType you entered was not unique. Fix the table!');
end

crystal = Crystals(I);
    