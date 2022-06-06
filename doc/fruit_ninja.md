# 水果忍者实验报告

软博21 沈冠霖 2021312593

## 1.问题转化

通过对偶化，这个问题可以转化为半空间求交问题。

对于每个竖直线段(a, b1),(a, b2)，其对偶化之后转化成两条平行线y=ax-b1,y=ax-b2中间的空隙。对于最终要求的直线y=ax-b，可以转化为点(a,b)。直线和线段相交，等价于存在线段上的一个点，使得直线对偶的点在这个点对偶得到的直线上，也就是直线对偶点位于线段对偶的两条平行线间的空隙中。这样，这个问题就转化成了n个竖直线段对应的2n条平行线划分的2n个半空间的求交问题。

## 2.算法实现

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