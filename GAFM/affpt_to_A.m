%
%  Copyright (c) 2018 James Pritts
%  Licensed under the MIT License (see LICENSE for details)
%
%  Written by James Pritts
%
function [A, A0] = affpt_to_A(affpt)
tmp = reshape([[affpt(:).x]; ...
               [affpt(:).y]; ...
               [affpt(:).a11]; ...
               [affpt(:).a12]; ...
               [affpt(:).a21]; ...
               [affpt(:).a22] ],6,[]);
geom = tmp([3 5 4 6 1 2],:);
A0 = reshape(geom,2,3,[]);
A = [ A0; ...
      repmat([0 0 1], [1 1 size(A0,3)])];
if ndims(size(A,3)) > 1
    A = squeeze(mat2cell(A,3,3,ones(1,size(A,3))));
end
if ndims(size(A0,3)) > 1
    A0 = squeeze(mat2cell(A0,2,3,ones(1,size(A0,3))));
end