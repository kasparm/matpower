function [A, l, u, vs, i1, iN] = retrieve_lin_constraints(om, name, idx)

if nargin < 3
    idx = {};
end
if isempty(idx)
    if numel(om.lin.idx.i1.(name)) == 1     %% simple named set
        A = om.lin.data.A.(name);
        l = om.lin.data.l.(name);
        u = om.lin.data.u.(name);
        if nargout > 3
            vs = om.lin.data.vs.(name);
            if nargout > 5
                i1 = om.lin.idx.i1.(name);      %% starting row index
                iN = om.lin.idx.iN.(name);      %% ending row index
            end
        end
    else                                    %% indexing required
        error('@opt_model/params_lin_constraint: linear constraint set ''%s'' requires an IDX arg', name);
    end
else
    % (calls to substruct() are relatively expensive ...
    % s = substruct('.', name, '{}', idx);
    % ... so replace it with these more efficient lines)
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
    A = subsref(om.lin.data.A, sc);
    l = subsref(om.lin.data.l, sc);
    u = subsref(om.lin.data.u, sc);
    if nargout > 3
        vs = subsref(om.lin.data.vs, sc);
        if nargout > 5
            sn = sc; sn(2).type = '()';         %% num array field
            i1 = subsref(om.lin.idx.i1, sn);    %% starting row index
            iN = subsref(om.lin.idx.iN, sn);    %% ending row index
        end
    end
end

end