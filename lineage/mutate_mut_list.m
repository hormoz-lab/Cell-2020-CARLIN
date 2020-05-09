function [relationship, mutations, type, preserved] = mutate_mut_list(parent_mut, child_mut)

    % Assume incompatible by default
    relationship = 'Incompatible';
    mutations = [];
    type = {};
    preserved = true(size(parent_mut));

    % Early special casing for template sequence
    if (isempty(parent_mut))
        if (isempty(child_mut))
            relationship = 'Identical';
            return;
        else
            relationship = 'Offspring';
            mutations = child_mut;
            type = repmat({'Additive'}, size(mutations));
            return;
        end
    else
        if (isempty(child_mut))
            return;
        end
    end
    
    N_mut_parent = length(parent_mut);
    N_mut_child  = length(child_mut);
    
    % Figure out which mutations are exclusive to parent/child, and which
    % are common
    parent_mut_exc = true(N_mut_parent,1);
    child_mut_exc  = true(N_mut_child, 1);        
    parent_idx=1;
    child_idx=1;
    
    while (parent_idx <= N_mut_parent && child_idx <= N_mut_child)
        if (isequal(parent_mut(parent_idx), child_mut(child_idx)))
            parent_mut_exc(parent_idx) = false;
            child_mut_exc(child_idx) = false;
            parent_idx = parent_idx+1;
            child_idx = child_idx+1;
        else
            if (parent_mut(parent_idx).loc_start < child_mut(child_idx).loc_start)
                parent_idx = parent_idx+1;
            elseif (parent_mut(parent_idx).loc_start > child_mut(child_idx).loc_start)
               child_idx = child_idx+1; 
            elseif (parent_mut(parent_idx).loc_start == child_mut(child_idx).loc_start)
                parent_idx = parent_idx+1;
                child_idx = child_idx+1;
            end
        end
    end
    
    preserved = ~parent_mut_exc;
    
    % If they don't have any exclusive mutations, they're identical
    if (~any(parent_mut_exc) && ~any(child_mut_exc))
        relationship = 'Identical';
        return;
    end
    
    % If the parent has any exclusive mutations, while the child has
    % none, they are incompatible (though not necessarily with reversed
    % arguments).
    if (any(parent_mut_exc) && ~any(child_mut_exc))
        return;
    end
    
    % If the child has any exclusive mutations, while the parent has none,
    % the child is a valid offspring, with purely additive mutations
    if (~any(parent_mut_exc) && any(child_mut_exc))
        mutations = child_mut(child_mut_exc);
        type = repmat({'Additive'}, size(mutations));
        relationship = 'Offspring';
        return;
    end
    
    % Comb through exclusive child mutations, and see if they can be placed
    % purely additively on the parent allele, or subtractively by subsuming
    % parental mutations.
    
    parent_mut_exc = find(parent_mut_exc);
    child_mut_exc  = find(child_mut_exc);
    
    N_mut_parent = length(parent_mut_exc);
    N_mut_child  = length( child_mut_exc);
    
    parent_idx=1;
    child_idx=1;    
    
    child_muts_to_apply = [];
    child_muts_to_apply_type = {};
    
    ref = CARLIN_def.getInstance;
    
    while (parent_idx <= N_mut_parent && child_idx <= N_mut_child)
        
        pm_start_loc_bp = parent_mut(parent_mut_exc(parent_idx)).loc_start;
        pm_end_loc_bp   = parent_mut(parent_mut_exc(parent_idx)).loc_end;
        cm_start_loc_bp = child_mut(child_mut_exc(child_idx)).loc_start;
        cm_end_loc_bp   = child_mut(child_mut_exc(child_idx)).loc_end;
        
        nested = (pm_start_loc_bp > cm_start_loc_bp && pm_end_loc_bp < cm_end_loc_bp);
        disjoint = (pm_start_loc_bp > cm_end_loc_bp);
                
        if (nested || disjoint)            
            pm_start_loc  = CARLIN_def.locate(ref, pm_start_loc_bp, ref.bounds.ordered);
            pm_end_loc    = CARLIN_def.locate(ref, pm_end_loc_bp  , ref.bounds.ordered);
            cm_start_loc  = CARLIN_def.locate(ref, cm_start_loc_bp, ref.bounds.ordered);
            cm_end_loc    = CARLIN_def.locate(ref, cm_end_loc_bp  , ref.bounds.ordered);
            parent_site_degraded = false;
            if (nested)
                if (pm_start_loc.rel == cm_start_loc.rel || pm_start_loc.rel==10 && strcmp(cm_start_loc.type, 'postfix'))
                    if ~((strcmp(pm_start_loc.type, 'pams') || strcmp(pm_start_loc.type, 'postfix')) && pm_start_loc.pos>3)
                        parent_site_degraded = true;
                    end
                end
                if (pm_end_loc.rel == cm_end_loc.rel || pm_end_loc.rel==10 && strcmp(cm_end_loc.type, 'postfix'))
                    if ~((strcmp(pm_end_loc.type, 'pams') || strcmp(pm_end_loc.type, 'postfix')) && pm_end_loc.pos>3)
                        parent_site_degraded = true;
                    end
                end
                if (parent_site_degraded)
                    return;
                end
                if (isempty(child_muts_to_apply) || child_muts_to_apply(end) ~= child_idx)
                    child_muts_to_apply      = [child_muts_to_apply; child_idx];
                    child_muts_to_apply_type = [child_muts_to_apply_type; {'Subtractive'}];
                end
                parent_idx = parent_idx+1;
            else
                if (pm_start_loc.rel == cm_end_loc.rel || pm_start_loc.rel==10 && strcmp(cm_end_loc.type, 'postfix'))
                    if ~((strcmp(pm_start_loc.type, 'pams') || strcmp(pm_start_loc.type, 'postfix')) && pm_start_loc.pos>3)
                        return;
                    end
                end
                if (isempty(child_muts_to_apply) || child_muts_to_apply(end) ~= child_idx)                
                    child_muts_to_apply = [child_muts_to_apply; child_idx];
                    child_muts_to_apply_type = [child_muts_to_apply_type; {'Additive'}];
                end
                child_idx = child_idx+1;
            end
        else
            return;
        end
    end
    
    % If there are remaining parent mutations, without any remaining child
    % mutations, incompatible.
    if (child_idx > N_mut_child)
        return;
    end
    
    if (child_idx < N_mut_child)    
        child_muts_to_apply      = [child_muts_to_apply; [child_idx+1:N_mut_child]'];
        child_muts_to_apply_type = [child_muts_to_apply_type; cell(N_mut_child-child_idx,1)];
        child_muts_to_apply_type(end-(N_mut_child-child_idx)+1:end) = {'Additive'}; 
    end

    type = child_muts_to_apply_type;
    mutations = child_mut(child_mut_exc(child_muts_to_apply));
    relationship = 'Offspring';
end