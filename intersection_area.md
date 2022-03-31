# 凸多边形求交实验报告

软博21 沈冠霖 2021312593

## 1.问题转化

凸多边形等价于其每条边构成的半空间的交集，根据凸集交运算的保凸性，不难得出，多个凸多边形的交集等价于这些凸多边形所有边构成的半空间的交集。因此，这个问题可以转化为求这些半空间的交集。

整体的算法我参考了这篇文章https://oi-wiki.org//geometry/half-plane/。其整体思路类似Graham Scan，算法如下：

1. 根据极角对所有线段进行排序，对于平行、共线的线段，只保留半空间最靠内的一个
2. 维护一个半空间-交点的双向队列，遍历排序好的线段。每进来一个新线段，先把队尾、队首所有在这个线段右侧的交点和半空间弹出，然后再和队尾的半空间求交点，加入队列
3. 最后用队首的线段去除队尾所有在其右侧的交点，然后再和队尾的半空间求交点，加入队列。这时候交点集合就是最终的交集了。
4. 使用向量叉乘公式求结果凸多边形的面积

极角排序的复杂度是O(nlogn)，队列处理中，由于除了队列首的每个线段只能进入队列1次，因此复杂度为O(n)，综合来看复杂度是O(nlogn)，满足题目数据范围的要求。

## 2.极角排序

Graham Scan的极角排序使用的是ToLeft方法，这里因为向量的极角范围是[0, 360)而不是Graham Scan的[0, 180)，因此不能直接用ToLeft进行求解。同时，虽然这里不需要考虑向量共线的问题，但是必须要处理好向量平行的问题。因此，我设计了如下的排序-筛选算法。

在排序阶段，我按照如下标准对向量进行比较：

1. 判断向量极角的范围。如果向量满足(y > 0 || (y == 0 && x > 0))，则其极角在[0, 180)范围内，否则在[180, 360)范围内。[0, 180)范围内的向量可以直接用ToLeft进行比较，而[180, 360)的向量将x， y都取反就转化成了[0, 180)的向量，可以ToLeft比较。使用范围比较也能够有效处理反向平行和反向共线的向量。
2. 对于在同一范围内的向量，使用ToLeft进行比较：a极角比b小当仅当ToLeft(o, a, b)严格成立（o是原点）。
3. 如果之前计算得到Area2=0，则说明两个向量同向平行或者同向共线，这个时候，我继续使用ToLeft去判断两者的关系：如果ToLeft(a.start, a.end, b.end)严格成立，则说明a和b同向平行，且b在内部，a应该被舍弃。如果这时候的ToLeft的Area2还是0，则说明a和b同向共线，那么任意舍弃一个即可。我们把要舍弃的排序在前面，要保留的排序在后面。

排序后需要做舍弃，我的方法是遍历排序后的向量，如果当前向量和其后一个向量同向平行或者同向共线，则说明当前向量需要舍弃，否则保留。判别同向平行或者同向共线，只需要判断向量极角在同一范围内(都在[0, 180)或[180, 360)范围内)，并且Area2(o, a, b)=0即可。

整体代码如下：

```C++
//排序的比较算法
bool operator<(const HalfSpace& b) const 
    {
        //judge up and down
        Point a_v = v;
        Point b_v = b.v;
        Point zero = Point(0.0, 0.0);
        bool whether_up_a = JudgeUp(a_v);
        bool whether_up_b = JudgeUp(b_v);
        if(whether_up_a == 0 && whether_up_b == 1)
        {
            return 0;
        }
        else if(whether_up_a == 1 && whether_up_b == 0)
        {
            return 1;
        }
        else if(whether_up_a == 0 && whether_up_b == 0)
        {
            a_v = zero - a_v;
            b_v = zero - b_v;
        }

        //judge to left
        int to_left = ToLeft(zero, a_v, b_v);
        if(to_left == 1)
        {
            return 1;
        }
        else if(to_left == -1)
        {
            return 0;
        }

        //handle parallel
        int to_left_points = ToLeft(s, t, b.t);
        if(to_left_points == 1)
        {
            return 1;
        }
        else 
        {
            return 0;
        }
    }

//同向平行/共线的判别算法
bool JudgeParallelColiniar(const HalfSpace& a, const HalfSpace& b)
{
    bool result = 0;
    bool a_up = JudgeUp(a.v);
    bool b_up = JudgeUp(b.v);
    if(a_up == b_up)
    {
        Point zero = Point(0.0, 0.0);
        double area = Area2(zero, a.v, b.v);
        if(area != 0)
        {
            result = 0;
        }
        else
        {
            result = 1;
        }
    }
    return result;
}

//整体代码
vector<HalfSpace> SortHalfSpaces(vector<HalfSpace>& half_spaces)
{
    vector<HalfSpace> sorted_list;
    sorted_list.clear();
    sort(half_spaces.begin(), half_spaces.end());
    for(int i = 0; i < half_spaces.size(); i ++)
    {
        if(i == half_spaces.size() - 1)
        {
            sorted_list.push_back(half_spaces[i]);
        }
        else 
        {
            HalfSpace current = half_spaces[i];
            HalfSpace next = half_spaces[i + 1];
            bool judge_result = !JudgeParallelColiniar(current, next);
            if(judge_result)
            {
                sorted_list.push_back(current);
            }
        }
    }
    return sorted_list;
}
```

## 3.主体流程

主体流程我参考了https://oi-wiki.org//geometry/half-plane/。与Graham算法相同的是，这里在进行极角排序后，由于向量有序，因此新进来的半空间只需要和队列首位比较即可。与Graham算法不同的是，因为Graham算法的极角都在[0, 180)范围内，因此只需要和队尾比较，而此问题极角范围变到了[0, 360)，**因此也需要和队首比较**。

1. 我分别维护半空间和交点两个双向队列。初始时候，半空间队列只有第一个半空间，交点队列为空。
2. 之后我遍历排序后的半空间列表，每加入一个新半空间，都要和交点-半空间队列的队首、队尾进行比较。如果被比较的交点在新的半空间**右侧**，说明这个交点已经**失效**了，需要出队列。
3. 排除所有失效交点后，我们再用新的半空间和队尾半空间求交点，两者分别加入两个队列。
4. 最终，因为最后一个加入队列的半空间只有一侧有交点，还要取出目前半空间队列首的半空间，在队尾进行失效交点检测--求交插入这个步骤。这样，交点队列构成的凸多边形就是半空间的交了。

整体代码如下：

```C++
vector<Point> HalfPlaneIntersection(vector<HalfSpace>& half_spaces)
{
    //initialize
    deque<Point> intersection_point_queue;
    deque<HalfSpace> half_space_queue; 
    intersection_point_queue.clear();
    half_space_queue.clear();
    vector<HalfSpace> sorted_half_spaces = SortHalfSpaces(half_spaces);
    half_space_queue.push_back(sorted_half_spaces[0]);

    //the main procedure
    for(int i = 1; i < sorted_half_spaces.size(); i ++)
    {
        HalfSpace current_half_space = sorted_half_spaces[i];
        //handle the back of the queue
        while(!intersection_point_queue.empty())
        {
            Point last_point = intersection_point_queue.back();
            int to_left = ToLeft(current_half_space.s, current_half_space.t, last_point);
            if(to_left < 0)
            {
                intersection_point_queue.pop_back();
                half_space_queue.pop_back();
            }
            else 
            {
                break;
            }
        }

        //handle the front of the queue
        while(!intersection_point_queue.empty())
        {
            Point first_point = intersection_point_queue.front();
            int to_left = ToLeft(current_half_space.s, current_half_space.t, first_point);
            if(to_left < 0)
            {
                intersection_point_queue.pop_front();
                half_space_queue.pop_front();
            }
            else 
            {
                break;
            }
        }
        HalfSpace last_half_space = half_space_queue.back();
        Point new_intersection_point = GetIntersectionPoint(current_half_space, last_half_space);
        half_space_queue.push_back(current_half_space);
        intersection_point_queue.push_back(new_intersection_point);
    }


    //use the queue front to handle the queue back
    HalfSpace current_half_space = half_space_queue.front();
    while(!intersection_point_queue.empty())
    {
        Point last_point = intersection_point_queue.back();
        int to_left = ToLeft(current_half_space.s, current_half_space.t, last_point);
        if(to_left < 0)
        {
            intersection_point_queue.pop_back();
            half_space_queue.pop_back();
         }
        else 
        {
            break;
        }
    }
    HalfSpace last_half_space = half_space_queue.back();
    Point new_intersection_point = GetIntersectionPoint(current_half_space, last_half_space);
    half_space_queue.push_back(current_half_space);
    intersection_point_queue.push_back(new_intersection_point);


    //push to the vector
    vector<Point> intersection_points;
    intersection_points.clear();
    while(!intersection_point_queue.empty())
    {
        Point first_point = intersection_point_queue.front();
        intersection_points.push_back(first_point);
        intersection_point_queue.pop_front();
    }
    return intersection_points;
}
```

## 4.凸多边形面积求解

这一部分我参考了https://blog.csdn.net/Mikchy/article/details/81490908。

凸多边形可以转化成若干个三角形的面积，而三角形面积可以用叉乘求解，因此可以用如下代码求解：

```C++
double ConvexPolygonArea(const vector<Point>& points)
{
    double result = 0.0;
    if(points.size() <= 2)
    {
        return 0.0;
    }
    Point start = points[0];
    for(int i = 1; i < points.size() - 1; i ++)
    {
        result += ((points[i] - start) ^ (points[i + 1] - start));
    }
    result = result / 2.0;
    return result;
}
```

## 5.踩坑

这个算法整体类似Graham Scan，但是也有不同。这里主要的坑点是加入新的半空间时，需要同时处理队首队尾，同时最后还要用队首再给队尾处理一次。需要这么做的原因是Graham Scan通过求LTL让极角范围变成了[0, 180)，而这里因为是向量极角排序，因此极角范围是[0, 360)，向量会绕一圈绕回来，因此需要首尾兼顾才可。