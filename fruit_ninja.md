# 水果忍者实验报告

软博21 沈冠霖 2021312593

## 1.问题转化

通过对偶化，这个问题可以转化为半空间求交问题。

对于每个竖直线段(a, b1),(a, b2)，其对偶化之后转化成两条平行线y=ax-b1,y=ax-b2中间的空隙。对于最终要求的直线y=ax-b，可以转化为点(a,b)。直线和线段相交，等价于存在线段上的一个点，使得直线对偶的点在这个点对偶得到的直线上，也就是直线对偶点位于线段对偶的两条平行线间的空隙中。这样，这个问题就转化成了n个竖直线段对应的2n条平行线划分的2n个半空间的求交问题。

## 2.算法实现

## 2.极角排序

Graham Scan的极角排序使用的是ToLeft方法，这里因为向量的极角范围是[0, 360)而不是Graham Scan的[0, 180)，因此不能直接用ToLeft进行求解。同时，虽然这里不需要考虑向量共线的问题，但是必须要处理好向量平行的问题。因此，我设计了如下的排序-筛选算法。

在排序阶段，我按照如下标准对向量进行比较：

1. 判断向量极角的范围。如果向量满足(y > 0 || (y == 0 && x > 0))，则其极角在[0, 180)范围内，否则在[180, 360)范围内。[0, 180)范围内的向量可以直接用ToLeft进行比较，而[180, 360)的向量将x， y都取反就转化成了[0, 180)的向量，可以ToLeft比较。使用范围比较也能够有效处理反向平行和反向共线的向量。

2. 对于在同一范围内的向量，使用ToLeft进行比较：a极角比b小当仅当ToLeft(o, a, b)严格成立（o是原点）。

   求。面积了。队列一个下：共线，这个时候，我继续使用ToLeft去判断两者的关系：如果ToLeft(a.start, a.end, b.end)严格成立，则说明a和b同向平行，且b在内部，a应该被舍弃。如果这时候的ToLeft的Area2还是0，则说明a和b同向共线，那么任意舍弃一个即可。我们把要舍弃的排序在前面，要保留的排序在后面。

半空间求交问题在1-4的PA中已经解决，我直接复用了那里的代码。需要修改的只有两点：

1.因为这里得到的半空间交可能不封闭，为了保证求交算法正确，我给整个空间添加了四个半空间作为边框（用四个横平竖直的半空间做约束，每个都保证能够覆盖题目给定的数据范围），具体代码如下：

```C++
        double smallest = -2000000000000.0;
        double biggest = 2000000000000.0;//比数据范围10^6的平方都大，必然能约束好整个空间
        Point down_left = Point(smallest, smallest);
        Point down_right = Point(biggest, smallest);
        Point up_right = Point(biggest, biggest);
        Point up_left = Point(smallest, biggest);
        HalfSpace l1 = HalfSpace(down_left, down_right);
        HalfSpace l2 = HalfSpace(down_right, up_right);
        HalfSpace l3 = HalfSpace(up_right, up_left);
        HalfSpace l4 = HalfSpace(up_left, down_left);
        half_spaces.push_back(l1);
        half_spaces.push_back(l2);
        half_spaces.push_back(l3);
        half_spaces.push_back(l4);
```

2.这里要求解的是半空间是否为空，按照我们的算法，如果最后得到的交点集合点数目>=3，则不为空，否则为空。