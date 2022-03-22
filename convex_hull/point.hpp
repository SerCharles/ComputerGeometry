#ifndef _POINT_
#define _POINT_
#include <vector>
#define MAX 2147483647
#define MIN -2147483648

class Point
{
public:
    int x;
    int y;
    Point() 
    {
        x = 0;
        y = 0;
    }
    Point(int a, int b)
    {
        x = a;
        y = b;
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
    q [2D Point]: [one of the three points]
    s [2D Point]: [one of the three points]
    
Returns:
    result [bool]: [1 if s is between p and q]
*/
bool Between(Point p, Point q, Point s)
{
    bool result = 0;
    int inner_product = (s.x - p.x) * (s.x - q.x) + (s.y - p.y) * (s.y - q.y);
    if(inner_product < 0)
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
    //TODO
    bool result = Area2(p, q, s);
    return result;
}

/*
Get the LTL(Lowest then Leftmost) point of a point list
Args:
    points [vector<Point>]: [the point list]
    
Returns:
    ltl [int]: [the id of the ltl point]
*/
int LTL(std::vector<Point> points)
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

#endif
