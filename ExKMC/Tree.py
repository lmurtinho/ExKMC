import numpy as np
import pandas as pd
from collections import Counter
from sklearn.cluster import KMeans
# removed . below to make it work on native Jupyter notebook
from splitters import get_min_mistakes_cut
from splitters import get_min_surrogate_cut

try:
    from graphviz import Source
    graphviz_available = True
except Exception:
    graphviz_available = False

BASE_TREE = ['IMM', 'NONE']

LEAF_DATA_KEY_X_DATA = 'X_DATA_KEY'
LEAF_DATA_KEY_Y = 'Y_KEY'
LEAF_DATA_KEY_X_CENTER_DOT = 'X_CENTER_DOT'
LEAF_DATA_KEY_SPLITTER = 'SPLITTER_KEY'


class Tree:

    def __init__(self, k, max_leaves=None, verbose=0, light=True,
                valid_col_idx=None, base_tree='IMM', n_jobs=None):
        """
        Constructor for explainable k-means tree.
        :param k: Number of clusters.
        :param max_leaves: Grow a tree with up to max_leaves. If None then a tree with k leaves will be constructed.
        :param verbose: Verbosity mode.
        :param light: If False, the object will store a copy of the input examples associated with each leaf.
        :pama valid_col_idx: array with indices of dimensions that can be used to split a tree. If None all dimensions can be used.
        :param base_tree: Specify weather the first k leaves are generated according to IMM splitting criteria or not. Valid values are ["IMM", "NONE"].
        :param n_jobs: The number of jobs to run in parallel.
        """
        self.k = k
        self.tree = None
        self.node_list = []
        self._leaves_data = {}
        self.max_leaves = k if max_leaves is None else max_leaves
        if self.max_leaves < k:
            raise Exception('max_trees must be greater or equal to k [%d < %d]' % (self.max_leaves, k))
        self.verbose = verbose
        self.light = light
        if base_tree not in BASE_TREE:
            raise Exception(base_tree + ' is not a supported base tree')
        self.base_tree = base_tree
        self.n_jobs = n_jobs if n_jobs is not None else 1
        self._feature_importance = None
        self.valid_col_idx = valid_col_idx

    def find_lca(self, x_data, y, valid_centers, dim):
        # print(valid_centers)
        centers = self.all_centers[valid_centers, dim]
        tree = self.trees_1d[dim]
        while True:
            if np.all(centers > tree.value):
                tree = tree.right
            elif np.all(centers <= tree.value):
                tree = tree.left
            else:
                return tree.value

    def find_largest_side(self, x_data, y, valid_centers, valid_cols):
        centers = self.all_centers[valid_centers]
        center_box = centers.max(axis=0) - centers.min(axis=0)
        val = center_box[valid_cols].max()
        return [i for i in range(len(center_box)) if center_box[i] == val][0]

    def find_fewest_mistakes(self, x_data, y, valid_centers, valid_cols):
        centers = self.all_centers[valid_centers]
        min_mistakes = np.inf
        best_dim = -1
        for dim in range(len(self.trees_1d)):
            if valid_cols[dim]:
                val = self.find_lca(x_data, y, valid_centers, dim)
                data_mask = x_data[:,dim] <= val
                assigns_mask = self.all_centers[y,dim] <= val
                mistakes_mask = data_mask != assigns_mask
                n_mistakes = mistakes_mask.sum()
                if n_mistakes < min_mistakes:
                    best_dim = dim
                    min_mistakes = n_mistakes
        return best_dim

    def find_mistakes_abv_avg(self, x_data, y, valid_centers, valid_cols):
        centers = self.all_centers[valid_centers]
        center_box = centers.max(axis=0) - centers.min(axis=0)
        mean_center_box = center_box.mean()
        below_mean_mask = center_box <= mean_center_box
        valid_cols[below_mean_mask] = False
        min_mistakes = np.inf
        best_dim = -1
        for dim in range(len(self.trees_1d)):
            if valid_cols[dim]:
                val = self.find_lca(x_data, y, valid_centers, dim)
                data_mask = x_data[:,dim] <= val
                assigns_mask = self.all_centers[y,dim] <= val
                mistakes_mask = data_mask != assigns_mask
                n_mistakes = mistakes_mask.sum()
                if n_mistakes < min_mistakes:
                    best_dim = dim
                    min_mistakes = n_mistakes
        return best_dim

    def build_from_1d(self, x_data, y, valid_centers,
                      find_feat="largest_side", find_val="lca"):
        # print(x_data.shape, y.shape, valid_centers.shape, valid_centers.sum())
        # print(self.all_centers[valid_centers])
        # print(valid_centers)

        node = Node()
        self.node_list.append(node)

        # base cases
        if x_data.shape[0] == 0:
            node.value = 0
            return node
        if valid_centers.sum() == 1:
            # print("single center")
            # print()
            node.value = np.argmax(valid_centers)
            return node
        if np.unique(y).shape[0] == 1:
            node.value = y[0]
            return node

        if find_feat == "largest_side":
            find_feat_func = self.find_largest_side
        elif find_feat == "fewest_mistakes":
            find_feat_func = self.find_fewest_mistakes
        elif find_feat == "fewest_mistakes_abv_avg":
            find_feat_func = self.find_mistakes_abv_avg
        if find_val == "lca":
            find_val_func = self.find_lca
        # check if valid_centers or self.all_centers should be passed
        # to functions
        valid_cols = np.array([1 if len(np.unique(x_data[:,i])) > 1 else 0
                                for i in range(x_data.shape[1])],
                                dtype=bool)
        node.feature = find_feat_func(x_data, y, valid_centers, valid_cols)
        node.value   = find_val_func(x_data, y, valid_centers, node.feature)
        # print(node.feature, node.value)
        # print()
        left_data_mask = x_data[:,node.feature] <= node.value
        left_assign_mask = self.all_centers[y,node.feature] <= node.value
        mistakes_mask = left_data_mask != left_assign_mask

        left_valid_centers = np.zeros(self.k, dtype=bool)
        right_valid_centers = np.zeros(self.k, dtype=bool)
        for i in range(self.k):
            if valid_centers[i]:
                if self.all_centers[i, node.feature] <= node.value:
                    left_valid_centers[i] = True
                else:
                    right_valid_centers[i] = True

        # left_valid_centers_mask = self.all_centers[valid_centers, node.feature] <= node.value
        # left_valid_centers = np.zeros(valid_centers.shape, dtype=bool)
        # left_valid_centers[valid_centers] = left_valid_centers_mask
        #
        # right_valid_centers = np.zeros(valid_centers.shape, dtype=bool)
        # right_valid_centers[valid_centers.astype(bool)] = ~left_valid_centers_mask

        node.left = self.build_from_1d(x_data[left_data_mask & ~mistakes_mask],
                                        y[left_data_mask & ~mistakes_mask],
                                        left_valid_centers,
                                        find_feat, find_val)
        node.right = self.build_from_1d(x_data[~left_data_mask & ~mistakes_mask],
                                         y[~left_data_mask & ~mistakes_mask],
                                         right_valid_centers,
                                         find_feat, find_val)

        return node


    def _build_tree(self, x_data, y, valid_centers, valid_cols):
        """
        Build a tree.
        :param x_data: The input samples.
        :param y: Clusters of the input samples, according to the kmeans classifier given (or trained) by fit method.
        :param valid_centers: Boolean array specifying which centers should be considered for the tree creation.
        :param valid_cols: Boolean array specifying which columns should be considered fot the tree creation.
        :return: The root of the created tree.
        """
        if self.verbose > 1:
            print('build node (samples=%d)' % x_data.shape[0])
        node = Node()
        self.node_list.append(node)
        if x_data.shape[0] == 0:
            node.value = 0
            return node
        elif valid_centers.sum() == 1:
            node.value = np.argmax(valid_centers)
            return node
        else:
            if np.unique(y).shape[0] == 1:
                node.value = y[0]
                return node
            else:

                # Verify data type is float64 prior to cython call
                x_data = x_data.astype(np.float64, copy=False)
                y = y.astype(np.int32, copy=False)
                self.all_centers = self.all_centers.astype(np.float64, copy=False)
                valid_centers = valid_centers.astype(np.int32, copy=False)
                valid_cols = valid_cols.astype(np.int32, copy=False)

                cut = get_min_mistakes_cut(x_data, y, self.all_centers, valid_centers, valid_cols, self.n_jobs)

                if cut is None:
                    node.value = np.argmax(valid_centers)
                else:
                    col = cut["col"]
                    threshold = cut["threshold"]
                    node.set_condition(col, threshold)

                    left_data_mask = x_data[:, col] <= threshold
                    matching_centers_mask = self.all_centers[:, col][y] <= threshold
                    mistakes_mask = left_data_mask != matching_centers_mask

                    left_valid_centers_mask = self.all_centers[valid_centers.astype(bool), col] <= threshold
                    left_valid_centers = np.zeros(valid_centers.shape, dtype=np.int32)
                    left_valid_centers[valid_centers.astype(bool)] = left_valid_centers_mask
                    right_valid_centers = np.zeros(valid_centers.shape, dtype=np.int32)
                    right_valid_centers[valid_centers.astype(bool)] = ~left_valid_centers_mask

                    node.left = self._build_tree(x_data[left_data_mask & ~mistakes_mask],
                                                 y[left_data_mask & ~mistakes_mask],
                                                 left_valid_centers,
                                                 valid_cols)
                    node.right = self._build_tree(x_data[~left_data_mask & ~mistakes_mask],
                                                  y[~left_data_mask & ~mistakes_mask],
                                                  right_valid_centers,
                                                  valid_cols)

                return node

    def build_1d_trees(self, x_data, y, type="cost_search"):

        if type == "cost_search":
            tree_func = self.tree_1d_cost_search
        elif type == "dasgupta":
            tree_func = self.tree_1d_cut_min
        self.trees_1d = [tree_func(x_data, y, i)
                         for i in range(x_data.shape[1])]

    def tree_1d_cut_min(self,x_data,y,dim):
        valid_cols = np.zeros(x_data.shape[1])
        valid_cols[dim] = 1
        return self._build_tree(x_data, y,
                                np.ones(self.all_centers.shape[0],
                                        dtype=np.int32),
                                valid_cols)

    def fit_from_1d(self,x_data,kmeans=None, type_1d="cost_search",
                    find_feat="largest_side", find_val="lca"):
        x_data = convert_input(x_data)

        if kmeans is None:
            if self.verbose > 0:
                print('Finding %d-means' % self.k)
            kmeans = KMeans(self.k, n_jobs=self.n_jobs, verbose=self.verbose, n_init=1, max_iter=40)
            kmeans.fit(x_data)
        else:
            assert kmeans.n_clusters == self.k

        y = np.array(kmeans.predict(x_data), dtype=np.int32)
        self.all_centers = np.array(kmeans.cluster_centers_, dtype=np.float64)
        valid_centers = np.ones(self.k, dtype=bool)

        self.build_1d_trees(x_data, y, type=type_1d)
        self.tree = self.build_from_1d(x_data, y, valid_centers,
                                        find_feat,find_val)

        self._feature_importance = np.zeros(x_data.shape[1])
        self.__fill_stats__(self.tree, x_data, y)

        return self

    def fit(self, x_data, kmeans=None):
        """
        Build a threshold tree from the training set x_data.
        :param x_data: The training input samples.
        :param kmeans: Trained model of k-means clustering over the training data.
        :return: Fitted threshold tree.
        """

        x_data = convert_input(x_data)

        if kmeans is None:
            if self.verbose > 0:
                print('Finding %d-means' % self.k)
            kmeans = KMeans(self.k, n_jobs=self.n_jobs, verbose=self.verbose, n_init=1, max_iter=40)
            kmeans.fit(x_data)
        else:
            assert kmeans.n_clusters == self.k

        y = np.array(kmeans.predict(x_data), dtype=np.int32)

        self.all_centers = np.array(kmeans.cluster_centers_, dtype=np.float64)

        if self.base_tree == "IMM":
            valid_cols = np.zeros(self.all_centers.shape[1], dtype=np.int32)
            valid_cols[self.valid_col_idx] = 1
            self.tree = self._build_tree(x_data, y,
                                         np.ones(self.all_centers.shape[0],
                                                 dtype=np.int32),
                                         valid_cols)
            leaves = self.k
        else:
            self.tree = Node()
            self.tree.value = 0
            leaves = 1

        if self.max_leaves > leaves:
            self.__gather_leaves_data__(self.tree, x_data, y)
            all_centers_norm_sqr = (np.linalg.norm(self.all_centers, axis=1) ** 2).astype(np.float64, copy=False)
            self.__expand_tree__(leaves, all_centers_norm_sqr)
            if self.light:
                self._leaves_data = {}

        self._feature_importance = np.zeros(x_data.shape[1])
        self.__fill_stats__(self.tree, x_data, y)

        return self

    def tree_1d_cost_search(self, x_data, y, dim):
        """
        Build a threshold tree from the training set x_data based on the
        algorithm to minimize the expected cost of searching through the tree,
        for a single dimension
        """
        idx_centers = self.all_centers[:,dim].argsort()
        centers_d   = self.all_centers[:,dim][idx_centers]
        assigns_d   = self.all_centers[y,dim]
        x_d         =  x_data[:,dim]

        separated = get_separations(x_d,centers_d,assigns_d)
        costs = separated.sum(axis=0)
        tree = tree_with_costs(costs, 0)
        tree_with_mids(tree,centers_d,idx_centers,dim)
        return tree

    def fit_predict(self, x_data, kmeans=None):
        """
        Build a threshold tree from the training set x_data, and returns the predicted clusters.
        :param x_data: The training input samples.
        :param kmeans: Trained model of k-means clustering over the training data.
        :return: The predicted clusters.
        """
        self.fit(x_data, kmeans)
        return self.predict(x_data)

    def predict(self, x_data):
        """
        Predict clusters for x_data.
        :param x_data: The input samples.
        :return: The predicted clusters.
        """
        x_data = convert_input(x_data)
        return self._predict_subtree(self.tree, x_data)

    def _predict_subtree(self, node, x_data):
        if node.is_leaf():
            return np.full(x_data.shape[0], node.value)
        else:
            ans = np.zeros(x_data.shape[0])
            left_mask = x_data[:, node.feature] <= node.value
            ans[left_mask] = self._predict_subtree(node.left, x_data[left_mask])
            ans[~left_mask] = self._predict_subtree(node.right, x_data[~left_mask])
            return ans

    def score(self, x_data):
        """
        Return the k-means cost of x_data.
        The k-means cost is the sum of squared distances of each point to the mean of points associated with the cluster.
        :param x_data: The input samples.
        :return: k-means cost of x_data.
        """
        x_data = convert_input(x_data)
        clusters = self.predict(x_data)
        cost = 0
        for c in range(self.k):
            cluster_data = x_data[clusters == c, :]
            if cluster_data.shape[0] > 0:
                center = cluster_data.mean(axis=0)
                cost += np.linalg.norm(cluster_data - center) ** 2
        return cost

    def surrogate_score(self, x_data):
        """
        Return the k-means surrogate cost of x_data.
        The k-means surrogate cost is the sum of squared distances of each point to the closest center of the kmeans given (or trained) in the fit method.
        k-means surrogate cost > k-means cost, as k-means cost is computed with respect to the optimal centers.
        :param x_data: The input samples.
        :return: k-means surrogate cost of x_data.
        """
        x_data = convert_input(x_data)
        clusters = self.predict(x_data)
        cost = 0
        for c in range(self.k):
            cluster_data = x_data[clusters == c, :]
            if cluster_data.shape[0] > 0:
                center = self.all_centers[c]
                cost += np.linalg.norm(cluster_data - center) ** 2
        return cost

    def _size(self):
        """
        Return the number of nodes in the threshold tree.
        :return: the number of nodes in the threshold tree.
        """
        return self.__size__(self.tree)

    def __size__(self, node):
        """
        Return the number of nodes in the subtree rooted by node.
        :param node: root of a subtree.
        :return: the number of nodes in the subtree rooted by node.
        """
        if node is None:
            return 0
        else:
            sl = self.__size__(node.left)
            sr = self.__size__(node.right)
            return 1 + sl + sr

    def _max_depth(self):
        """
        Return the depth of the threshold tree.
        :return: the depth of the threshold tree.
        """
        return self.__max_depth__(self.tree)

    def __max_depth__(self, node):
        """
        Return the depth of the subtree rooted by node.
        :param node: root of a subtree.
        :return: the depth of the subtree rooted by node.
        """
        if node is None:
            return -1
        else:
            dl = self.__max_depth__(node.left)
            dr = self.__max_depth__(node.right)
            return 1 + max(dl, dr)

    def plot(self, filename="test", feature_names=None):
        if not graphviz_available:
            raise Exception("Required package is missing. Please install graphviz")

        if self.tree is not None:
            dot_str = ["digraph ClusteringTree {\n"]
            queue = [self.tree]
            nodes = []
            edges = []
            id = 0
            while len(queue) > 0:
                curr = queue.pop(0)
                if curr.is_leaf():
                    label = "%s\nsamples=\%d\nmistakes=\%d" % (str(curr.value), curr.samples, curr.mistakes)
                else:
                    feature_name = curr.feature if feature_names is None else feature_names[curr.feature]
                    condition = "%s <= %.3f" % (feature_name, curr.value)
                    label = "%s\nsamples=\%d" % (condition, curr.samples)
                    queue.append(curr.left)
                    queue.append(curr.right)
                    edges.append((id, id + len(queue) - 1))
                    edges.append((id, id + len(queue)))
                nodes.append({"id": id,
                              "label": label,
                              "node": curr})
                id += 1
            for node in nodes:
                dot_str.append("n_%d [label=\"%s\"];\n" % (node["id"], node["label"]))
            for edge in edges:
                dot_str.append("n_%d -> n_%d;\n" % (edge[0], edge[1]))
            dot_str.append("}")
            dot_str = "".join(dot_str)
            try:
                s = Source(dot_str, filename=filename + '.gv', format="png")
                s.view()
            except Exception as e:
                print(dot_str)
                raise e

    def __expand_tree__(self, size, all_centers_norm_sqr):
        if size < self.max_leaves:
            if self.verbose > 1:
                print('expand tree. size %d/%d' % (size, self.max_leaves))

            best_splitter = None
            leaf_to_split = None
            leaf_count = 1
            for leaf in self._leaves_data:
                if self.verbose > 1:
                    print('-- expand leaf. %d/%d (samples=%d)' % (
                        leaf_count, len(self._leaves_data), self._leaves_data[leaf][LEAF_DATA_KEY_X_DATA].shape[0]))
                if LEAF_DATA_KEY_SPLITTER not in self._leaves_data[leaf]:
                    self._leaves_data[leaf][LEAF_DATA_KEY_SPLITTER] = self.__expand_leaf__(leaf, all_centers_norm_sqr)
                leaf_splitter = self._leaves_data[leaf][LEAF_DATA_KEY_SPLITTER]
                if leaf_splitter is not None:
                    if best_splitter is None:
                        best_splitter = leaf_splitter
                        leaf_to_split = leaf
                    elif leaf_splitter["cost_gain"] < best_splitter["cost_gain"]:
                        best_splitter = leaf_splitter
                        leaf_to_split = leaf
                leaf_count += 1
            if best_splitter is not None:
                col = best_splitter["col"]
                threshold = best_splitter["threshold"]
                self.__split_leaf__(leaf_to_split,
                                    col,
                                    threshold,
                                    best_splitter["center_left"],
                                    best_splitter["center_right"])

                X = self._leaves_data[leaf_to_split][LEAF_DATA_KEY_X_DATA]
                y = self._leaves_data[leaf_to_split][LEAF_DATA_KEY_Y]
                X_center_dot = self._leaves_data[leaf_to_split][LEAF_DATA_KEY_X_CENTER_DOT]
                left_mask = X[:, col] <= threshold

                del self._leaves_data[leaf_to_split]

                self._leaves_data[leaf_to_split.left] = {LEAF_DATA_KEY_X_DATA: X[left_mask],
                                                         LEAF_DATA_KEY_Y: y[left_mask],
                                                         LEAF_DATA_KEY_X_CENTER_DOT: X_center_dot[left_mask]}
                self._leaves_data[leaf_to_split.right] = {LEAF_DATA_KEY_X_DATA: X[~left_mask],
                                                          LEAF_DATA_KEY_Y: y[~left_mask],
                                                          LEAF_DATA_KEY_X_CENTER_DOT: X_center_dot[~left_mask]}
                self.__expand_tree__(size + 1, all_centers_norm_sqr)

    def __gather_leaves_data__(self, node, x_data, y):
        if node.is_leaf():
            self._leaves_data[node] = {LEAF_DATA_KEY_X_DATA: x_data,
                                       LEAF_DATA_KEY_Y: y,
                                       LEAF_DATA_KEY_X_CENTER_DOT: np.dot(x_data, self.all_centers.T).astype(np.float64,
                                                                                                             copy=False)}
        else:
            left_mask = x_data[:, node.feature] <= node.value
            self.__gather_leaves_data__(node.left, x_data[left_mask], y[left_mask])
            self.__gather_leaves_data__(node.right, x_data[~left_mask], y[~left_mask])

    def __expand_leaf__(self, leaf, all_centers_norm_sqr):
        leaf_data = self._leaves_data[leaf]
        mistakes_counter = Counter([curr_y for curr_y in leaf_data[LEAF_DATA_KEY_Y] if curr_y != leaf.value])
        if len(mistakes_counter) == 0:
            return None

        # Verify data type is float64 prior to cython call
        X = leaf_data[LEAF_DATA_KEY_X_DATA].astype(np.float64, copy=False)
        X_center_dot = leaf_data[LEAF_DATA_KEY_X_CENTER_DOT].astype(np.float64, copy=False)
        all_centers_norm_sqr = all_centers_norm_sqr.astype(np.float64, copy=False)

        min_cut = get_min_surrogate_cut(X, X_center_dot, X_center_dot.sum(axis=0), all_centers_norm_sqr, self.n_jobs)

        if min_cut is not None:
            pre_split_cost = self.__get_leaf_pre_split_cost__(X_center_dot, all_centers_norm_sqr)
            splitter = {"col": min_cut["col"],
                        "threshold": min_cut["threshold"],
                        "cost_gain": min_cut["cost"] - pre_split_cost,
                        "center_left": min_cut["center_left"],
                        "center_right": min_cut["center_right"]}
            return splitter
        else:
            return None

    def __get_leaf_pre_split_cost__(self, X_center_dot, all_centers_norm_sqr):
        n = X_center_dot.shape[0]

        cost_per_center = (n * all_centers_norm_sqr) - 2 * X_center_dot.sum(axis=0)
        best_center = cost_per_center.argmin()
        return cost_per_center[best_center]

    def __split_leaf__(self, leaf, feature, value, left_cluster, right_cluster):
        leaf.feature = feature
        leaf.value = value

        leaf.left = Node()
        leaf.left.value = left_cluster

        leaf.right = Node()
        leaf.right.value = right_cluster

    def __fill_stats__(self, node, x_data, y):
        node.samples = x_data.shape[0]
        if not node.is_leaf():
            self._feature_importance[node.feature] += 1
            left_mask = x_data[:, node.feature] <= node.value
            self.__fill_stats__(node.left, x_data[left_mask], y[left_mask])
            self.__fill_stats__(node.right, x_data[~left_mask], y[~left_mask])
        else:
            node.mistakes = len([cluster for cluster in y if cluster != node.value])

    def feature_importance(self):
        return self._feature_importance


class Node:
    def __init__(self):
        self.feature = None
        self.value = None
        self.samples = None
        self.mistakes = None
        self.left = None
        self.right = None

    def is_leaf(self):
        return (self.left is None) and (self.right is None)

    def set_condition(self, feature, value):
        self.feature = feature
        self.value = value


def convert_input(data):
    if isinstance(data, list):
        data = np.array(data, dtype=np.float64)
    elif isinstance(data, np.ndarray):
        data = data.astype(np.float64, copy=False)
    elif isinstance(data, pd.DataFrame):
        data = data.values.astype(np.float64, copy=False)
    else:
        raise Exception(type(data) + ' is not supported type')
    return data

def find_first_above(val, seq):
    if not len(seq):
        return 0
    elif len(seq) == 1:
        return int(val >= seq[0])
    mid = len(seq) // 2
    if seq[mid] > val:
        return find_first_above(val, seq[:mid])
    else:
        return mid + 1 + find_first_above(val, seq[mid+1:])


def get_separations(data, centers, assigns):
    mids = (centers[1:] + centers[:-1]) / 2
    # print(centers, mids)
    above_d = np.array([find_first_above(d, mids) for d in data])
    above_c_dict = {c: find_first_above(c, mids) for c in centers}
    above_c = np.array([above_c_dict[a] for a in assigns])
    return np.array([[int( (i >= min(j)) and (i < max(j)) )
                         for i in range(len(mids))]
                        for j in zip(above_d, above_c)])

def tree_with_costs(costs, prev):
    if not len(costs):
        return Node()
    min_cost = costs.min()
    if not min_cost:
        cut = costs.argmin()
        # print("clean cut identified", cut)
        node = Node()
        node.feature = cut + prev
        # node.mistakes = 0
        node.left  = tree_with_costs(costs[:cut], prev)
        node.right = tree_with_costs(costs[cut+1:], node.feature+1)
        return node
    costs = costs / min_cost
    n = len(costs)
    j=1
    lc=rc=0
    t = np.ceil(np.log2(n)**(1/3))
    r = max(np.ceil(1 - 1/t), 2)
    # print("n", n, "t", t, "r", r)
#     print(costs, j, lc, rc, r, t)
    return tree_with_pos_costs(costs,j,lc,rc,r,t,prev)

def tree_with_pos_costs(costs,j,lc,rc,r,t,prev):
    # print("infos:", costs, j, lc, rc, r, t)
    if not len(costs):
        # print()
        return Node()
    min_cost = r**(j-1)
    max_cost = r**j
    mask = (costs >= min_cost) * (costs <= max_cost)
#     print("Step 1")
    while not np.any(mask):
#         print(min_cost, max_cost)
        j += 1
        min_cost = r**(j-1)
        max_cost = r**j
        mask = (costs >= min_cost) * (costs <= max_cost)
        lc = rc = 0 # TODO: isso fica na segunda versão?
    idx_mask = np.where(mask)[0]
    # print("idx_mask:", idx_mask)

    add = 1
    if lc == t:
        lc = 0
        x_idx = idx_mask[0]
    elif rc == t:
        rc = 0
        x_idx = idx_mask[-1]
    else:
        add = 0
        m_start = idx_mask[0]
        m_end = idx_mask[-1]
        M = costs[m_start:m_end+1]
        # print("M:", M)
        x_idx = (m_start + m_end) // 2
        step = 0
        while not mask[x_idx]:
            step += 1
            x_idx += step
            if mask[x_idx]:
                break
            step += 1
            x_idx -= step

    node = Node()
    node.feature = x_idx + prev
    # node.mistakes = costs[x_idx]
    # print("result:", node.feature, node.value)
    # print()
    # TODO: check valores de lc e rc na segunda versão
    node.left  = tree_with_pos_costs(costs[:x_idx],1,lc+add,0,r,t,prev)
    node.right = tree_with_pos_costs(costs[x_idx+1:],1,0,rc+add,r,t,node.feature+1)
    return node

def tree_with_mids(tree,centers,indices,dim):
    mids = (centers[1:] + centers[:-1]) / 2
    nodes = [(None, None, tree)]
    while nodes:
        cut, side, node = nodes.pop(0)
        if node.left:
            nodes.append((node.feature, "left", node.left))
            nodes.append((node.feature, "right", node.right))
            node.value = mids[node.feature]
            node.feature = dim
        elif side == "left":
            node.value = indices[cut]
        else:
            node.value = indices[cut+1]
