%
%  Copyright (c) 2018 James Pritts
%  Licensed under the MIT License (see LICENSE for details)
%
%  Written by James Pritts
%
function u = affpt_to_pt3x3(affpt)
u = A_to_pt3x3(affpt_to_A(affpt));
