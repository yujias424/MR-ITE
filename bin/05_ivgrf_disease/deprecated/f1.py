def convert_color(color):
    try:
        color = pl.get_cmap(color)
    except:
        pass
    
    if color == "shap_red":
        color = colors.red_rgb
    elif color == "shap_blue":
        color = colors.blue_rgb
    
    return color

def convert_ordering(ordering, shap_values):
    if issubclass(type(ordering), OpChain):
        ordering = ordering.apply(Explanation(shap_values))
    if issubclass(type(ordering), Explanation):
        if "argsort" in [op["name"] for op in ordering.op_history]:
            ordering = ordering.values
        else:
            ordering = ordering.argsort.flip.values
    return ordering

def get_sort_order(dist, clust_order, cluster_threshold, feature_order):
    """ Returns a sorted order of the values where we respect the clustering order when dist[i,j] < cluster_threshold
    """
    
    clust_inds = np.argsort(clust_order)

    feature_order = feature_order.copy()
    for i in range(len(feature_order)-1):
        ind1 = feature_order[i]
        next_ind = feature_order[i+1]
        next_ind_pos = i + 1
        for j in range(i+1,len(feature_order)):
            ind2 = feature_order[j]

            if dist[ind1,ind2] <= cluster_threshold:
                
                if dist[ind1,next_ind] > cluster_threshold or clust_inds[ind2] < clust_inds[next_ind]:
                    next_ind = ind2
                    next_ind_pos = j
    
        for j in range(next_ind_pos, i+1, -1):
            feature_order[j] = feature_order[j-1]
        feature_order[i+1] = next_ind
    
    return feature_order

def merge_nodes(values, partition_tree):
    """ This merges the two clustered leaf nodes with the smallest total value.
    """
    M = partition_tree.shape[0] + 1

    ptind = 0
    min_val = np.inf
    for i in range(partition_tree.shape[0]):
        ind1 = int(partition_tree[i,0])
        ind2 = int(partition_tree[i,1])
        if ind1 < M and ind2 < M:
            val = np.abs(values[ind1]) + np.abs(values[ind2])
            if val < min_val:
                min_val = val
                ptind = i
                #print("ptind", ptind, min_val)

    ind1 = int(partition_tree[ptind,0])
    ind2 = int(partition_tree[ptind,1])
    if ind1 > ind2:
        tmp = ind1
        ind1 = ind2
        ind2 = tmp
    
    partition_tree_new = partition_tree.copy()
    for i in range(partition_tree_new.shape[0]):
        i0 = int(partition_tree_new[i,0])
        i1 = int(partition_tree_new[i,1])
        if i0 == ind2:
            partition_tree_new[i,0] = ind1
        elif i0 > ind2:
            partition_tree_new[i,0] -= 1
            if i0 == ptind + M:
                partition_tree_new[i,0] = ind1
            elif i0 > ptind + M:
                partition_tree_new[i,0] -= 1

            
        if i1 == ind2:
            partition_tree_new[i,1] = ind1
        elif i1 > ind2:
            partition_tree_new[i,1] -= 1
            if i1 == ptind + M:
                partition_tree_new[i,1] = ind1
            elif i1 > ptind + M:
                partition_tree_new[i,1] -= 1
    partition_tree_new = np.delete(partition_tree_new, ptind, axis=0)

    # update the counts to be correct
    fill_counts(partition_tree_new)
    
    return partition_tree_new, ind1, ind2

def fill_counts(partition_tree):
    """ This updates the 
    """
    M = partition_tree.shape[0] + 1
    for i in range(partition_tree.shape[0]):
        val = 0
        if partition_tree[i,0] < M:
            ind = int(partition_tree[i,0])
            val += 1
        else:
            ind = int(partition_tree[i,0])-M
            val += partition_tree[ind,3]
        if partition_tree[i,1] < M:
            ind = int(partition_tree[i,1])
            val += 1
        else:
            ind = int(partition_tree[i,1])-M
            val += partition_tree[ind,3]
        partition_tree[i,3] = val

def sort_inds(partition_tree, leaf_values, pos=None, inds=None):
    if inds is None:
        inds = []
    
    if pos is None:
        partition_tree = fill_internal_max_values(partition_tree, leaf_values)
        pos = partition_tree.shape[0]-1
    
    M = partition_tree.shape[0] + 1
        
    if pos < 0:
        inds.append(pos + M)
        return
    
    left = int(partition_tree[pos, 0]) - M
    right = int(partition_tree[pos, 1]) - M
    
    
    left_val = partition_tree[left,3] if left >= 0 else leaf_values[left + M]
    right_val = partition_tree[right,3] if right >= 0 else leaf_values[right + M]
    
    if left_val < right_val:
        tmp = right
        right = left
        left = tmp
    
    sort_inds(partition_tree, leaf_values, left, inds)
    sort_inds(partition_tree, leaf_values, right, inds)
    
    return inds(base)