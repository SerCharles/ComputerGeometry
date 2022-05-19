#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stack>
#include <deque>
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
};

/*
Get the signed area of a triangle defined by three points

Args:
    p [2D Point]: [one of the three points of the triangle]
    q [2D Point]: [one of the three points of the triangle]
    s [2D Point]: [one of the three points of the triangle]
    
Returns:
    area [double]: [the signed area * 2 of the triangle, it is positive if s is at the left of the line pq]
*/
double Area2(const Point& p, const Point& q, const Point& s)
{
    double area = p.x * q.y - p.y * q.x + q.x * s.y - q.y * s.x + s.x * p.y - s.y * p.x;
    return area;
}

/*
The toLeft test, testing the relative position of point s and line pq

Args:
    p [2D Point]: [one of the three points of the triangle]
    q [2D Point]: [one of the three points of the triangle]
    s [2D Point]: [one of the three points of the triangle]
    
Returns:
    result [int]: [1 if s is at the left of line pq, -1 if s is at the right of line pq, 0 if the three points are coliniar]
*/
int ToLeft(const Point& p, const Point& q, const Point& s)
{
    double area2 = Area2(p, q, s);
    int result = 0;
    if(area2 > 0)
    {
        result = 1;
    }
    else if(area2 < 0)
    {
        result = -1;
    }
    else 
    {
        result = 0;
    }
    return result;
}

/*
Judge whether a point's polar angle is in[0, 180), iff (y > 0 || (y == 0 && x > 0))

Args:
    a [2D Point]: [the point]

Returns:
    result [bool]: [whether the point's polar angle is in[0, 180), 1 if so]
*/
bool JudgeUp(const Point& a)
{
    bool result = 0;
    if (a.y > 0 || (a.y == 0 && a.x > 0))
    {
        result = 1;
    }
    return result;
}


/*
Definition of a half space
Args:
    s [double]: the start point of the segment of the half space
    t [double]: the end point of the segment of the half space
    v [double]: the vector of the half space, t - s
*/
class HalfSpace
{
public:
    Point s;
    Point t;
    Point v;
    HalfSpace() {}
    HalfSpace(Point a, Point b)
    {
        s = a;
        t = b;
        v = t - s;
    }    

    /*
        Judge the order of two half spaces in the sorting algorithm
        The criterions: 
            1.For polar angles in [0, 360), half spaces with smaller polar angle should be ahead
            2.If two half spaces have the same polar angle, leave the inner one behind.
            3.If two half spaces are coliniar and of the same direction, everything is ok since they are the same
        The algorithm of judging order:
            1.Judge the range of the polar angles, whether they are in [0, 180) or not
            2.Judge toleft of (0, v1, v2), consider their polar angles, if in [180, 360), v1 := -v1, v2 := -v2
            3.If they are parallel, judge toleft of (s1, t1, t2)
        The algorithm of sifting:
            for any line in points[i], if lines[i + 1] is not parallel and not coliniar, then lines[i] should be saved
    */
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
};

/*
    If the two halfspaces are of the same direction, judge whether they are parallel or coliniar
    for any line in points[i], if lines[i + 1] is not parallel and not coliniar, then lines[i] should be saved

    Args:
        a [2D HalfSpace]: [one of the half spaces]
        b [2D HalfSpace]: [one of the half spaces]

    Returns:
        result [bool]: [1 if they are parallel or coliniar, 0 if else]
*/
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

/*
Get the Intersection point of two lines

Args:
    a [2D HalfSpace]: [line a]
    b [2D HalfSpace]: [line b]

Returns:
    result [2D Point]: [the intersection point]
*/
Point GetIntersectionPoint(const HalfSpace& a, const HalfSpace& b)
{
    double c1 = (b.t - b.s) ^ (a.s - b.s);
    double c2 = (b.t - b.s) ^ (a.t - b.s);
    double new_x = (a.s.x * c2 - a.t.x * c1) / (c2 - c1);
    double new_y = (a.s.y * c2 - a.t.y * c1) / (c2 - c1);
    Point result = Point(new_x, new_y);
    return result;
}

/*
Get the area of a convex polygon, which is defined by a list of points in counter clockwise direction

Args:
    points [vector of Point]: [the points defining the polygon]

Returns:
    result [double]: [the area]
*/
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

/*
Sort the half spaces by polar angle, and remove redundant half spaces
The criterions: 
    1.For polar angles in [0, 360), half spaces with smaller polar angle should be ahead
    2.If two half spaces have the same polar angle, leave the inner one behind.
    3.If two half spaces are coliniar and of the same direction, everything is ok since they are the same
The algorithm of judging order:
    1.Judge the range of the polar angles, whether they are in [0, 180) or not
    2.Judge toleft of (0, v1, v2), consider their polar angles, if in [180, 360), v1 := -v1, v2 := -v2
    3.If they are parallel, judge toleft of (s1, t1, t2)
The algorithm of sifting:
    for any line in points[i], if lines[i + 1] is not parallel and not coliniar, then lines[i] should be saved

Args:
    half_spaces [vector of HalfSpace]: [the original half space list]

Returns:
    sorted_list [vector of HalfSpace]: [the sorted unique half space list]
*/
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

/*
The main algorithm of Half Plane Intersection

Args:
    half_spaces [vector of HalfSpace]: [the original half space list]

Returns:
    intersection_points [vector of Point]: [the intersection points forming the intersection area in counter clockwise direction]
*/
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

int main()
{
    //input
    vector<HalfSpace> half_spaces;
    half_spaces.clear();
    int polygons = 0;
    scanf("%d", &polygons);
    for(int i = 0; i < polygons; i ++)
    {
        int num = 0;
        scanf("%d", &num);
        vector<Point> points;
        points.clear();
        for(int j = 0; j < num; j ++)
        {
            int x = 0;
            int y = 0;
            scanf("%d %d", &x, &y);
            Point new_point = Point(double(x), double(y));
            points.push_back(new_point);
        }
        for(int j = 0; j < points.size(); j ++)
        {
            int next = (j + 1) % points.size();
            HalfSpace new_halfspace = HalfSpace(points[j], points[next]);
            half_spaces.push_back(new_halfspace);
        }
    }
    
    //output
    vector<Point> result_points = HalfPlaneIntersection(half_spaces);
    double result = ConvexPolygonArea(result_points);
    printf("%.4lf", result);
    return 0;
}