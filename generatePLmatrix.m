function plMatrix = generatePLmatrix(chMatrix)
% generatePLmatrix function
%
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% Generate path loss matrix
% Input: chMatrix
%
% Output: plMatrix

plMatrix = nan(size(chMatrix));

for snapshotIdx = 1 : size(chMatrix,3)
    for TxIdx = 1 : size(chMatrix,1)
        for RxIdx = 1 : size(chMatrix,2)

            if size(chMatrix{TxIdx,RxIdx,snapshotIdx},1) > 0
                plMatrix(TxIdx,RxIdx,snapshotIdx) = -mag2db(sum(abs(chMatrix{TxIdx,RxIdx,snapshotIdx}.h)));
            end


        end
    end
end

end