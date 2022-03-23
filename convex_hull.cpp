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
    x [int]: the x coordinate of the point
    y [int]: the y coordinate of the point
    id [int]: the id of the point, starting from 1 to n
*/
class Point
{
public:
    int x;
    int y;
    int id;
    Point() 
    {
        x = 0;
        y = 0;
        id = 0;
    }
    Point(int a, int b, int i)
    {
        x = a;
        y = b;
        id = i;
    }

};

/*
Get the signed area of a triangle defined by three points

Args:
    p [2D Point]: [one of the three points of the triangle]
    q [2D Point]: [one of the three points of the triangle]
    s [2D Point]: [one of the three points of the triangle]
    
Returns:
    area [int]: [the signed area * 2 of the triangle, it is positive if s is at the left of the line pq]
*/
int Area2(Point p, Point q, Point s)
{
    int area = p.x * q.y - p.y * q.x + q.x * s.y - q.y * s.x + s.x * p.y - s.y * p.x;
    return area;
}

/*
The between test, testing whether s is between p and q, if the three points are coliniar

Args:
    p [2D Point]: [one of the three points]
    s [2D Point]: [one of the three points]
    q [2D Point]: [one of the three points]
    
Returns:
    result [bool]: [1 if s is between p and q]
*/
bool Between(Point p, Point s, Point q)
{
    bool result = 0;
    int inner_product = (p.x - s.x) * (s.x - q.x) + (p.y - s.y) * (s.y - q.y);
    if(inner_product > 0)
    {
        result = 1;
    }
    return result;
}

/*
The between test, testing whether p, q, s are coliniar

Args:
    p [2D Point]: [one of the three points]
    s [2D Point]: [one of the three points]
    q [2D Point]: [one of the three points]
    
Returns:
    result [bool]: [1 if the three points are coliniar)
*/
bool Coliniar(Point p, Point s, Point q)
{
    bool result = 0;
    int area2 = Area2(p, q, s);
    if(area2 == 0)
    {
        result = 1;
    }
    return result; 
}

/*
The toLeft test, testing the relative position of point s and line pq

Args:
    p [2D Point]: [one of the three points of the triangle]
    q [2D Point]: [one of the three points of the triangle]
    s [2D Point]: [one of the three points of the triangle]
    
Returns:
    result [bool]: [1 if s is at the left of line pq]
*/
bool ToLeft(Point p, Point q, Point s)
{
    int area2 = Area2(p, q, s);
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
int LTL(std::vector<Point>& points)
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
The main function of Graham Scan algorithm to get the convex hull
Args:
    points [vector<Point>]: [the point list]
    
Returns:
    S [vector<Point>]: [the polar points of the convex hull]
*/
vector<Point> GrahamScan(vector<Point>& points)
{
    //sort the points by polar angle
    int ltl = LTL(points);
    Swap(points, 0, ltl);
    PolarAngle polar_angle_comparator(points[0]);
    sort(points.begin() + 1, points.end(), polar_angle_comparator);

    /*
    printf("---------------------\n");
    for(int i = 0; i < points.size(); i ++)
    {
        printf("%d %d %d\n", points[i].x, points[i].y, points[i].id);
    }*/

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
            S.push_back(w);
            T.pop_back();
        }
        else 
        {
            S.pop_back();
        }
    }

    //additional effort: the extreme points in the last edge is not identified
    vector<bool> visit;
    visit.clear();
    for(int i = 0; i <= points.size(); i ++)
    {
        visit.push_back(0);
    }
    for(int i = 0; i < S.size(); i ++)
    {
        int id = S[i].id;
        visit[id] = 1;
    }
    Point first = S[0];
    Point last = S[S.size() - 1];
    for(int i = 0; i < points.size(); i ++)
    {
        int id = points[i].id;
        if(visit[id] == 0)
        {
            bool coliniar = Coliniar(first, last, points[i]);
            if(coliniar)
            {
                S.push_back(points[i]);
            }
        }
    }

    /*
    printf("---------------------\n");
    for(int i = 0; i < S.size(); i ++)
    {
        printf("%d %d %d\n", S[i].x, S[i].y, S[i].id);
    }
    */
    return S;
}

int main()
{
    //input
    int total_number = 0;
    int big_number = 1000000007;
    vector<Point> total_points;
    total_points.clear();
    scanf("%d", &total_number);
    for(int i = 1; i <= total_number; i ++)
    {
        int x = 0, y = 0;
        scanf("%d %d", &x, &y);
        Point new_point(x, y, i);
        total_points.push_back(new_point);
    }

    //output
    vector<Point> polar_points = GrahamScan(total_points);
    long long int result = 1;
    for(int i = 0; i < polar_points.size(); i ++)
    {
        result = (result * polar_points[i].id) % big_number;
    }
    result = (result * polar_points.size()) % big_number;
    printf("%lld", result);
    return 0;
}