//整体算法参考了论文Computing the Convex Hull of Line Intersections 
//感谢陈嘉杰同学提供论文资料
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stack>
using namespace std;


/*
Definition of a point
Args:
    x [double]: the x coordinate of the point
    y [double]: the y coordinate of the point
*/
class Point
{
public:
    double x;
    double y;
    Point() 
    {
        x = 0;
        y = 0;
    }
    Point(double a, double b)
    {
        x = a;
        y = b;
    }
    Point operator+(const Point& a) const
	{
        double new_x = x + a.x;
        double new_y = y + a.y;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
    Point operator-(const Point& a) const
	{
        double new_x = x - a.x;
        double new_y = y - a.y;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	void operator+=(const Point& a)
	{
		x += a.x;
        y += a.y;
	}
	void operator-=(const Point& a)
	{
		x -= a.x;
        y -= a.y;
	}
	Point operator*(const double& a) const
	{
		double new_x = x * a;
        double new_y = y * a;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	Point operator/(const double& a) const
	{
		double new_x = x / a;
        double new_y = y / a;
        Point new_point = Point(new_x, new_y);
		return new_point;
	}
	void operator*=(const double& a)
	{
		x *= a;
        y *= a;	
    }
	void operator/=(const double& a)
	{
        x /= a;
        y /= a;	
    }
    double operator^(const Point& b) const 
    {
        double result = x * b.y - y * b.x;
        return result;
    }
    //sort by LTL, used in removing redundant
    bool operator<(const Point& b) const 
    {
        bool result = 0;
        if(y < b.y)
        {
            result = 1;
        }
        else if(y == b.y)
        {
            if(x < b.x)
            {
                result = 1;
            }
        }
        return result;
    }

    bool operator==(const Point& b) const 
    {
        double dist = GetDistance(b);
        if(dist <= 1e-12)
        {
            return 1;
        }
        else 
        {
            return 0;
        }
    }

    double GetDistance(const Point& b) const 
    {
        double dist = sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y));
        return dist;
    }
    
};

class Line
{
public:
    double k;
    double b;
    Line()
    {
        k = 0;
        b = 0;
    }
    Line(double kk, double bb)
    {
        k = kk;
        b = bb;
    }


    //Sort the Lines by k, used in getting the points
    bool operator<(const Line& line) const 
    {
        bool result = 0;
        if (k < line.k)
        {
            result = 1;
        }
        else if(k == line.k)
        {
            if(b < line.b)
            {
                result = 1;
            }
        }
        return result;
    }

    /*
    Get the intersection between two lines
    Args:
        line [Line]: [the other line, which is not parallel to this one]
        
    Returns:
        result [point]: [the intersecting point]
    */
    Point GetIntersection(const Line& line) const 
    {
        double x = -(b - line.b) / (k - line.k);
        double y = k * x + b;
        if(x == -0)
        {
            x = 0.0;
        }
        if(y == -0)
        {
            y = 0.0;
        }
        Point result = Point(x, y);
        return result;
    }
};

/*
Get the signed area of a triangle defined by three points

Args:
    p [Point]: [one of the three points of the triangle]
    q [Point]: [one of the three points of the triangle]
    s [Point]: [one of the three points of the triangle]
    
Returns:
    area [double]: [the signed area * 2 of the triangle, it is positive if s is at the left of the line pq]
*/
double Area2(const Point& p, const Point& q, const Point& s)
{
    double area = p.x * q.y - p.y * q.x + q.x * s.y - q.y * s.x + s.x * p.y - s.y * p.x;
    return area;
}


/*
The between test, testing whether s is between p and q, if the three points are coliniar

Args:
    p [Point]: [one of the three points]
    s [Point]: [one of the three points]
    q [Point]: [one of the three points]
    
Returns:
    result [bool]: [1 if s is between p and q]
*/
bool Between(Point p, Point s, Point q)
{
    bool result = 0;
    double inner_product = (p.x - s.x) * (s.x - q.x) + (p.y - s.y) * (s.y - q.y);
    if(inner_product > 0)
    {
        result = 1;
    }
    return result;
}

/*
The coliniar test, testing whether p, q, s are coliniar

Args:
    p [Point]: [one of the three points]
    s [Point]: [one of the three points]
    q [Point]: [one of the three points]
    
Returns:
    result [bool]: [1 if the three points are coliniar)
*/
bool Coliniar(Point p, Point s, Point q)
{
    bool result = 0;
    double area2 = Area2(p, q, s);
    if(area2 == 0)
    {
        result = 1;
    }
    return result; 
}

/*
The toLeft test, testing the relative position of point s and line pq

Args:
    p [Point]: [one of the three points of the triangle]
    q [Point]: [one of the three points of the triangle]
    s [Point]: [one of the three points of the triangle]
    
Returns:
    result [bool]: [1 if s is at the left of line pq]
*/
bool ToLeft(Point p, Point q, Point s)
{
    double area2 = Area2(p, q, s);
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

/*
Get the LTL(Lowest then Leftmost) point of a point list
Args:
    points [vector<Point>]: [the point list]
    
Returns:
    ltl [int]: [the index of the ltl point in the array]
*/
int LTL(vector<Point>& points)
{
    int ltl = 0;
    for(int i = 1; i < points.size(); i ++)
    {
        if(points[i].y < points[ltl].y || (points[i].y == points[ltl].y && points[i].x < points[ltl].x))
        {
            ltl = i;
        }
    }
    return ltl;
}

/*
The criterion used in sorting the points by polarangle
criterion: a < b iff b is at the left of the line base-a
Args:
    base [Point]: the LTL point of the points
*/
struct PolarAngle
{
    Point base;
    PolarAngle() {}
    PolarAngle(Point b)
    {
        base = b;
    }
    bool operator() (Point a, Point b) 
    { 
        bool result = ToLeft(base, a, b);
        return result;
    }
};

/*
Swap the ith and jth points
Args:
    points [vector<Point>]: [the point list]
    i [int]: [the index of the first point to be swapped]
    j [int]: [the index of the second point to be swapped]
*/
void Swap(vector<Point>& points, int i, int j)
{
    Point temp;
    temp = points[i];
    points[i] = points[j];
    points[j] = temp;
}

/*
Get the useful points for constructing the convex hull of lines
Args:
    lines [vector<Line>]: [the input lines]
Returns:
    points [vector<Point>]: [the useful points for constructing the convex hull of lines]
*/
vector<Point> GetUsefulPoints(vector<Line>& lines)
{
    vector<Point> points;
    points.clear();
    sort(lines.begin(), lines.end());
    for(int i = 0; i < lines.size(); i ++)
    {
        Line line = lines[i];
        double k = line.k;
        double next_k = k;
        int j = (i + 1) % lines.size();
        while(j != i)
        {
            double new_k = lines[j].k;
            bool get_intersection_point = 0;
            bool end_loop = 0;

            //judge status
            //parallel to original: pass
            if(new_k == k)
            {
                get_intersection_point = 0;
                end_loop = 0;
            }
            else 
            {
                //a new intersecting line
                if(next_k == k)
                {
                    next_k = new_k;
                    get_intersection_point = 1;
                    end_loop = 0;                    
                }
                else 
                {   
                    //a new intersecting line, but parallel to the previous one
                    if(new_k == next_k)
                    {
                        get_intersection_point = 1;
                        end_loop = 0;      
                    }
                    else 
                    {
                        get_intersection_point = 0;
                        end_loop = 1;
                    }
                }
            }

            if(get_intersection_point)
            {
                Point new_point = line.GetIntersection(lines[j]);
                points.push_back(new_point);    
            }
            if(end_loop)
            {
                break;
            }
            else 
            {
                j = (j + 1) % lines.size();
            }
        }
    }
    return points;
}

/*
Remove the redundant points of a convex hull

Args:
    raw_points [vector<Point>]: [the points]

Returns:
    points [vector<Point>]: [the points without redundant]
*/
vector<Point> RemoveRedundantPoints(vector<Point>& raw_points)
{
    sort(raw_points.begin(), raw_points.end());
    vector<bool> save_list;
    vector<Point> points;
    save_list.clear();
    points.clear();
    int n = raw_points.size();
    for(int i = 0; i < n; i ++)
    {
        save_list.push_back(1);
    }
    for(int i = 0; i < n - 1; i ++)
    {
        if(raw_points[i] == raw_points[i + 1])
        {
            save_list[i + 1] = 0;
        }
    }
    for(int i = 0; i < n; i ++)
    {
        if(save_list[i])
        {
            points.push_back(raw_points[i]);
        }
    }
    return points;
}

/*
Graham Scan algorithm to get the convex hull
Args:
    points [vector<Point>]: [the point list]
Returns:
    S [vector<Point>]: [the result point list]
*/
vector<Point> GrahamScan(vector<Point>& points)
{
    //sort the points by polar angle
    int ltl = LTL(points);
    Swap(points, 0, ltl);
    PolarAngle polar_angle_comparator(points[0]);
    sort(points.begin() + 1, points.end(), polar_angle_comparator);

    //initialize the stacks
    //stack top: the last one of the vector
    //push_back = push stack, pop_back = pop stack
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
   return S;
}


/*
Get the edge points of a convex hull

Args:
    points [vector<Point>]: [the points on the convex hull, at least 3 points]

Returns:
    edge_points [vector<Point>]: [the edge points on the convex hull]
*/
vector<Point> GetEdgePoints(vector<Point> &points)
{
    int n = points.size();
    vector<bool> save_list;
    vector<Point> edge_points;
    save_list.clear();
    edge_points.clear();
    for(int i = 0; i < n; i ++)
    {
        save_list.push_back(1);
    }
    for(int i = 0; i < n; i ++)
    {
        int current = i;
        int next = (i + 1) % n;
        int last = (i - 1 + n) % n;
        bool coliniar = Coliniar(points[last], points[current], points[next]);
        if(coliniar)
        {
            save_list[current] = 0;
        }
    }
    for(int i = 0; i < n; i ++)
    {
        if(save_list[i])
        {
            edge_points.push_back(points[i]);
        }
    }
    return edge_points;
}

/*
Get the max distance of two points by brute force

Args:
    points [vector<Point>]: [the points on the convex hull]

Returns:
    max_dist [double]: [the max distance of two points]
*/
double GetMaxDistance(vector<Point>& points)
{
    double max_dist = 0.0;
    int n = points.size();
    for(int i = 0; i < n; i ++)
    {
        for(int j = 0; j < n; j ++)
        {
            double dist = points[i].GetDistance(points[j]);
            if(dist > max_dist)
            {
                max_dist = dist;
            }
        }
    }
    return max_dist;
}

/*
Judge if all lines are parallel, use this to avoid degeneration
Args:
    lines [vector<Line>]: [the lines]
Return:
    result [bool]: [if all parallel, return 1, else return 0]
*/
double JudgeAllParallel(vector<Line>& lines)
{
    bool result = 1;
    for(int i = 0; i < lines.size() - 1; i ++)
    {
        if(lines[i + 1].k != lines[i].k)
        {
            result = 0;
            break;
        }
    }
    return result;
}


int main()
{
    vector<Line> lines;
    lines.clear();
    int n;
    scanf("%d", &n);
    for(int i = 0; i < n; i ++)
    {
        int k, b;
        scanf("%d%d", &k, &b);
        Line new_line = Line(double(k), double(b));
        lines.push_back(new_line);
    }
    if(JudgeAllParallel(lines))
    {
        printf("%lf", 0.0);
        return 0;
    }

    double max_dist;
    vector<Point> raw_points = GetUsefulPoints(lines);
    vector<Point> points = RemoveRedundantPoints(raw_points);
    if(points.size() <= 1)
    {
        max_dist = 0.0;
    }
    else if(points.size() == 2)
    {
        max_dist = points[0].GetDistance(points[1]);
    }
    else 
    {
        vector<Point> convex_hull_points = GrahamScan(points);
        vector<Point> edge_points = GetEdgePoints(convex_hull_points);
        max_dist = GetMaxDistance(edge_points);
    }
    printf("%lf", max_dist);
    return 0;
}

