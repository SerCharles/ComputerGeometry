# 凸包实验报告

软博21 沈冠霖 2021312593

## 1.算法

我使用了Graham Scan算法，对于任意情况，它都能在O(nlogn)时间下完成凸包计算，对于10^5个点，这种速度足够了。

Graham Scan算法的流程如下：首先找到LTL(Lowest then Leftmost)点，以这个点为基准对所有其他点按照极角进行排序。根据凸集分离定理，LTL点必为极点（而且必然是极边的端点）。排序的方法是进行ToLeft test，p在q前当且仅当ToLeft(LTL, p, q)==True。

之后用两个栈，S，T进行操作，具体代码如下，流程和课件相同：

```C++
	vector<Point> T;
    vector<Point> S;
    T.clear();
    S.clear();
    S.push_back(points[0]);
    S.push_back(points[1]);
    for(int i = points.size() - 1; i >= 2; i --)
    {
        T.push_back(points[i]);
    }

    //main procedure of Graham Scan
    while(!T.empty())
    {
        int t_size = T.size();
        int s_size = S.size();
        Point w = T[t_size - 1];
        Point v = S[s_size - 1];
        Point u = S[s_size - 2];
        bool to_left = ToLeft(u, v, w);
        if(to_left)
        {
            Point t_top = T[T.size() - 1];
            S.push_back(t_top);
            T.pop_back();
        }
        else 
        {
            S.pop_back();
        }
    }
```

对于ToLeft，我采用了课件的实现，能有效解决多点共线问题和水平、竖直线问题

```c++
long long Area2(Point p, Point q, Point s)
{
    long long area = p.x * q.y - p.y * q.x + q.x * s.y - q.y * s.x + s.x * p.y - s.y * p.x;
    return area;
}
bool Between(Point p, Point s, Point q)
{
    bool result = 0;
    long long inner_product = (p.x - s.x) * (s.x - q.x) + (p.y - s.y) * (s.y - q.y);
    if(inner_product > 0)
    {
        result = 1;
    }
    return result;
}
bool ToLeft(Point p, Point q, Point s)
{
    long long area2 = Area2(p, q, s);
    bool result = 0;
    if(area2 > 0)
    {
        result = 1;
    }
    else if(area2 < 0)
    {
        result = 0;
    }
    else 
    {
        result = Between(p, q, s);
    }
    return result;
}
```

这种实现方法能够让极角排序阶段的点集构成偏序关系（任意两点都能比较出大小），把最后一条极边以外的极点都正确提取。

但是这种方法，对于最后一条极边，只能提取出两个端点，中间的点无法提取，因为最后一条极边我们希望排序后得到的顺序是由远及近，而这种排序算法会由近及远（其他边都是由近及远）。因此我们最后还要进行特判，遍历所有的点，保留位于最后一条极边两个端点连线上的点。因为本题不需要按顺序输出，所以仅需要一次遍历，O(n)的复杂度。如果需要按顺序输出，还需要进行一次排序，变成O(nlogn)复杂度，都不影响整体复杂度。

## 2.踩坑

首先是最后一条极边的问题，这个需要额外考虑。

其次是long long的问题，除了最后求解需要long long之外，前面的x y坐标也需要long long。。。因为最大是10^5，所以用int的话Area2会溢出。感谢张雄帅同学的提醒，这个花费了我四五个小时的时间。。。