# 范围查询实验报告

软博21 沈冠霖 2021312593

## 1.算法介绍

我使用的是课件中提到的Range Tree算法完成这道题目。Range Tree算法是构造两段BBST（平衡二叉搜索树），先在X方向搜索点，再在Y方向搜索点。这个算法分为两部分：一部分是构造Range Tree数据结构，另一部分是给定范围进行查询。

## 2.Range Tree的构造

### 2.1 算法流程

先是X方向的BBST构造。为了实现方便，我定义的X方向BBST与课件略有不同，每个BBST Node有如下属性：

- 是否是叶子结点leaf，如果当前Node的点个数小于等于50，就是叶子结点
- 点列表points
- 点个数point_number
- 代表值x，代表这个Node内所有点的x坐标中位数
- 左子树left_son，其内部所有点的x坐标都小于等于代表值x
- 右子树right_son，其内部所有点的x坐标都大于等于代表值x
- Y方向的BBST树y_tree

X方向定义代码如下，Y方向与之类似不再赘述：

```c++
class RangeNodeX
{
public:
    int depth = 0; //the depth of current node
    int point_number = 0; //the number of points in the range
    int x = 0; //the x value of the node 
    bool leaf = 0; //whether it is leaf node or not
    vector<Point> points; //the points in the range 
    RangeNodeY* y_tree = NULL; //the y tree
    RangeNodeX* left_son = NULL;
    RangeNodeX* right_son = NULL;

    ~RangeNodeX()
    {
        this->points.clear();
        delete(this->y_tree);
        delete(this->left_son);
        delete(this->right_son);
    }
    RangeNodeX(int depth_, vector<Point>& points_);
};
```

我使用递归的方式构造X方向的BBST，流程如下：

1. 存储输入的点，记录点个数等。
2. 使用C++ STL库自带的nth_element()函数求顶点列表的中位数，参考了http://c.biancheng.net/view/7476.html。假设有k个点，这个函数可以让第k/2个点的值变成中位数的值，并且让其前面的点的值都小于等于它，后面的点的值都大于等于它。其原理是反复调用快排中的pivot点扫描，虽然其最坏复杂度是O(n^2)，但是平均而言，其期望复杂度为O(n)，符合我们的要求。
3. 判断是否是叶子节点，如果是就停止，否则按照中位数来切分顶点列表，构建左右子树。
4. 所有子树构建完毕之后，构造Y方向的BBST树。

X方向构造代码如下，Y方向与之类似不再赘述：

```c++
RangeNodeX::RangeNodeX(int depth_, vector<Point>& points_)
{
    //set values
    this->depth = depth_;
    this->point_number = points_.size();
    this->points.clear();
    for(int i = 0; i < this->point_number; i ++)
    {
        this->points.push_back(points_[i]);
    }

    //get median and split the list
    int median_index = this->point_number / 2;
    //参考了https://www.cnblogs.com/zzzlight/p/14298223.html#
    nth_element(this->points.begin(), this->points.begin() + median_index, this->points.end(), CompareX);
    this->x = this->points[median_index].x;

    //leaf node 
    if(this->point_number <= SPLIT_THRESHOLD)
    {
        this->leaf = 1;
        this->left_son = NULL;
        this->right_son = NULL;
        this->y_tree = NULL;
    }
    //build son
    else 
    {
        this->leaf = 0;
        vector<Point> left_list;
        vector<Point> right_list;
        left_list.clear();
        right_list.clear();
        for(int i = 0; i < median_index; i ++)
        {
            left_list.push_back(this->points[i]);
        }

        for(int i = median_index; i < this->point_number; i ++)
        {
            right_list.push_back(this->points[i]);
        }
        this->left_son = new RangeNodeX(this->depth + 1, left_list);
        this->right_son = new RangeNodeX(this->depth + 1, right_list);
        left_list.clear();
        right_list.clear();

        //build Y
        this->y_tree = new RangeNodeY(1, this, 0, this->point_number);
    }  
}

```

之后是Y方向的BBST构造。每个BBST Node有如下属性：

- 是否是叶子结点leaf，如果当前节点的点个数小于等于50，就是叶子结点
- 对应的X方向BBST树father
- 在X方向点列表的起始位置start
- 在X方向点列表占用的点个数point_number
- 代表值y，代表这个Node内所有点的y坐标中位数
- 左子树left_son，其内部所有点的y坐标都小于等于代表值y
- 右子树right_son，其内部所有点的y坐标都大于等于代表值y

我使用递归的方式构造X方向的BBST，流程如下：

1. 存储输入的起始位置、点个数等。
2. 用上文提到的nth_element函数求其对应X方向点列表在[start, start + point_number)间的中位数，同时把这个范围的点修改顺序，使其满足左边的Y<=中位数Y<=右边的Y这一性质。
3. 判断是否是叶子节点，如果是就停止，否则按照中位数来切分顶点列表，构建左右子树。这里就不需要实际去复制内存之类的了，只需要算出左右子树的起始位置和顶点个数即可。

### 2.2 复杂度分析

首先是分析其层数：因为求中位数能保证左右子树点差距不超过1，可以认为X方向、Y方向的树层数都是O(logn)。

之后是分析N个节点情况下，Y子树的时空复杂度。根节点构造Y子树的主定理公式为T(n)=2T(n/2)+O(n)，时间复杂度是O(nlogn)。Y子树每个节点的空间复杂度是O(1)，主定理公式为T(n)=2T(n/2)+O(1)，因此空间复杂度是O(n)。

之后是分析构造X子树的时空复杂度。X子树构造分为四步，第一步存储消耗时间复杂度O(n)，空间复杂度O(n)；第二步求中位数消耗时间复杂度O(n)，空间复杂度O(1)；第三步切分左右子树消耗时间复杂度O(n)，空间复杂度O(n)；第四步构造Y子树消耗时间复杂度O(nlogn)，空间复杂度O(n)。其时间主定理公式为T(n)=2T(n/2)+O(nlogn)，时间复杂度是O(n*(logn)^2)。空间主定理公式为T(n)=2T(n/2)+O(n)，空间复杂度是O(nlogn)。

综上所述，构造Range Tree的时间复杂度是O(n*(logn)^2)，空间复杂度是O(nlogn)。这种方法构造的效率并不是很高，没有达到课件提出的O(nlogn)的时间复杂度，但是对于本题来说够用了。同时，为了提高效率，我没有让节点个数达到1才变成叶子结点，我设置了节点个数小于等于50都是叶子结点，这样可以有效降低层数，提高实际的运行效率。

## 3.查询

对于每个矩形包围盒，我的查询流程如下：

1. 首先在X子树下查询，寻找矩形定义的下限x_min和上限x_max的搜索路径route_left和route_right。由于层数是O(logn)，因此时空复杂度都是O(logn)。
2. 其次求这两条路径的最近公共祖先（LCA)，时空复杂度都是O(logn)。
3. 之后是求正则子集（canonical subset ），从LCA开始遍历两条路径。如果左路径上的Node，其x值大于等于x_min，则说明其右子节点的所有x值大于等于x_min，将其右子节点加入。如果右路径上的Node,其x值小于等于x_max，则说明其左子节点的所有x值小于等于x_max，将其左子节点加入。对于叶子结点，直接加入。时空复杂度都是O(logn)。
4. 最后是对于每个正则子集对应的Node，在Y子树下查询。在Y子树下查询的1、2、3步骤与X子树下完全一致，也是先找搜索路径，再求最近公共祖先，再求正则子节点，不过这里求完正则子节点就可以直接加和了。时空复杂度都是O(logn)。共有O(logn)个正则子集对应的Node。

综上所述，查询的时空复杂度都是O((logn)^2)，对于m次查询，其时空复杂度就是O(m*(logn)^2)。

X方向查询代码如下，Y方向与之类似不再赘述：

```c++
int CountPoints(BoundingBox& box, RangeNodeX* root)
{
    int min_x = box.left;
    int max_x = box.right;
    vector<RangeNodeX*> route_left;
    vector<RangeNodeX*> route_right;
    vector<RangeNodeX*> canonical_nodes;
    route_left.clear();
    route_right.clear();
    canonical_nodes.clear();

    //get the left route
    RangeNodeX* current_left = root;
    route_left.push_back(current_left);
    while(current_left->leaf == 0)
    {
        if(min_x <= current_left->x)
        {
            current_left = current_left->left_son;
            route_left.push_back(current_left);
        }
        else 
        {
            current_left = current_left->right_son;
            route_left.push_back(current_left);
        }
    }

    //get the right route
    RangeNodeX* current_right = root;
    route_right.push_back(current_right);
    while(current_right->leaf == 0)
    {
        if(max_x >= current_right->x)
        {
            current_right = current_right->right_son;
            route_right.push_back(current_right);
        }
        else 
        {
            current_right = current_right->left_son;
            route_right.push_back(current_right);
        }
    }

    //get lca
    int lca = 0;
    for(int i = 0; i < min(route_left.size(), route_right.size()); i ++)
    {
        if(route_left[i] == route_right[i])
        {
            lca = i;
        }
        else 
        {
            break;
        }
    }

    //get canonical nodes
    if(route_left[lca]->leaf)
    {
        canonical_nodes.push_back(route_left[lca]);
    }
    else 
    {
        for(int i = lca + 1; i < route_left.size(); i ++)
        {
            if(route_left[i]->leaf)
            {
                canonical_nodes.push_back(route_left[i]);
            }
            else if(route_left[i]->x >= min_x)
            {
                canonical_nodes.push_back(route_left[i]->right_son);
            }
        }
        for(int i = lca + 1; i < route_right.size(); i ++)
        {
            if(route_right[i]->leaf)
            {
                canonical_nodes.push_back(route_right[i]);
            }
            else if(route_right[i]->x <= max_x)
            {
                canonical_nodes.push_back(route_right[i]->left_son);
            }
        }
    }
    
    //count
    int result = 0;
    for(int i = 0; i < canonical_nodes.size(); i ++)
    {
        int the_result = 0;
        if(canonical_nodes[i]->leaf)
        {
            for(int j = 0; j < canonical_nodes[i]->point_number; j ++)
            {
                if(box.JudgeInside(canonical_nodes[i]->points[j]))
                {
                    the_result += 1;
                }
            }
        }
        else 
        {
            the_result = CountPoints(box, canonical_nodes[i]->y_tree);
        }
        result = result + the_result;
    }
    return result;
}
```

## 4.总结与反思

整体算法的时间复杂度是O((m+n)* * (logn)^2)，空间复杂度是O(m*(logn)^2 + nlogn)。n是150000，m是100000，log(n)是12左右，对于这道题，使用这种算法要更佳，因为如果使用KD-Tree，其时间复杂度是O(m* * sqrt(n))，sqrt(n)是400左右，比(logn)^2要大一些。

我的整体实现也并不是很完美：一方面，复杂度还是偏高，使用fractional cascading 能把单次查询时间复杂度降到O(logn)，构造的时空复杂度降到O(nlogn)。另一方面，由于X和Y方向的构造、搜索流程大致相同，由于时间原因，我没来得及做很好的封装，导致重复代码较多，代码难以比较好的维护。想要构造能长期使用的底层库，需要付出更多的时间精力去优化算法和代码结构。