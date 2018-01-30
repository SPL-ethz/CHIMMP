%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ETH Zurich, Switzerland
% UC Santa Barbara, USA
%
% Project:  Novartis (UCSB+ETHZ)
% Year:     2018
% MATLAB:   R2017a, Windows 10
% Authors:  Dave Ochsenbein (DRO)
%
% Purpose:
% This function serves to replace the Matlab toolbox function normc. It
% should behave identically, though performance might be somewhat worse
% (not tested).
%
% Last modified:
% - 30. Jan 2018, DRO: Initial creation
%
%
% Input arguments:
% - M:                      (m x n)-matrix (unnormalized)
%
% Output argument:
% - M_normed:               (m x n)-matrix (normalized)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2018 David Ochsenbein
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

function [ M_normed ] = normc( M )

    M_normed = sqrt(M.^2 ./ sum(M.^2)) .* sign(M);

end

