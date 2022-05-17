#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;


//DCEL data structure
class Vertex;
class HalfEdge;
class Face;

class Point
{
public:
    int id;
    long long x;
    long long y;
    Face* bucket = NULL; //the bucket of the point
    Point(int id_, long long x_, long long y_)
    {
        this->id = id;
        this->x = x_;
        this->y = y_;
    }
    ~Point(){}
};

class Vertex
{
public:
    Point* p; //the coordinate of the vertex
    HalfEdge* incident_halfedge; //reference to the first outgoing incident half-edge
    Vertex(Point* p_)
    {
        this->p = p_;
        this->incident_halfedge = NULL;
    }
    ~Vertex(){}

};

class HalfEdge 
{
public:
    HalfEdge* twin = NULL; //reference to the twin half edge
    Vertex* origin = NULL; //reference to the origin vertex
    Face* incident_face = NULL; //reference to the left incident face
    HalfEdge* previous_halfedge = NULL; //reference to CCW previous half-edge
    HalfEdge* successor_halfedge = NULL; //reference to CCW next half-edge
    HalfEdge(){}
    ~HalfEdge(){}

};

class Face
{
public:
    HalfEdge* incident_halfedge = NULL; //the first half edge incident to the face from left
    vector<Point*> bucket_points; //the unassigned points of the bucket
    Face()
    {
        this->bucket_points.clear();
    }
    ~Face()
    {
        this->bucket_points.clear();
    }
}

//useful functions
/*
Get the signed area of a triangle defined by three points

Args:
    p [Point*]: [one of the three points of the triangle]
    q [Point*]: [one of the three points of the triangle]
    s [Point*]: [one of the three points of the triangle]
    
Returns:
    area [long long]: [the signed area * 2 of the triangle, it is positive if s is at the left of the line pq]
*/
long long Area2(Point* p, Point* q, Point* s)
{
    long long area = p->x * q->y - p->y * q->x + q->x * s->y - q->y * s->x + s->x * p->y - s->y * p->x;
    return area;
}

/*
The between test, testing whether s is between p and q, if the three points are coliniar

Args:
    p [Point*]: [one of the three points]
    s [Point*]: [one of the three points]
    q [Point*]: [one of the three points]
    
Returns:
    result [bool]: [1 if s is between p and q]
*/
bool Between(Point* p, Point* s, Point* q)
{
    bool result = 0;
    long long inner_product = (p->x - s->x) * (s->x - q->x) + (p->y - s->) * (s->y - q->y);
    if(inner_product > 0)
    {
        result = 1;
    }
    return result;
}

/*
The coliniar test, testing whether p, q, s are coliniar

Args:
    p [Point*]: [one of the three points]
    s [Point*]: [one of the three points]
    q [Point*]: [one of the three points]
    
Returns:
    result [bool]: [1 if the three points are coliniar)
*/
bool Coliniar(Point* p, Point* s, Point* q)
{
    bool result = 0;
    long long area2 = Area2(p, q, s);
    if(area2 == 0)
    {
        result = 1;
    }
    return result; 
}

/*
The toLeft test, testing the relative position of point s and line pq

Args:
    p [Point*]: [one of the three points of the triangle]
    q [Point*]: [one of the three points of the triangle]
    s [Point*]: [one of the three points of the triangle]
    
Returns:
    result [bool]: [1 if s is at the left of line pq]
*/
bool ToLeft(Point* p, Point* q, Point* s)
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

//the main class
class DelaunayTriangulation
{
    vector<Point*> points; //all the points
    vector<Point*> unassigned_points; //the unassigned points
    vector<Face*> faces; //the faces of the entire place, including the immense one
    vector<Vertex*> vertexs; //all vertexs, used in recording and deletion
    vector<Vertex*> half_edges; //all half edges, used in recording and deletion

    DelaunayTriangulation(vector<Point>& points_)
    {
        //record all the points
        this->points.clear();
        this->unassigned_points.clear();
        this->faces.clear();
        this->vertexs.clear();
        this->half_edges.clear();
        for(int i = 0; i < points_.size(); i ++)
        {
            Point* new_point = new Point(points_[i].id, points_[i].x, points_[i].y);
            this->points.push_back(new_point);
        }

        //set the first triangle
        //use the first and second point as point p and q, find a non coliniar point as point r
        Point* p = points[0];
        Point* q = points[1];
        Point* r = NULL;
        int r_id = -1;
        for(int i = 2; i < this->points.size(); i ++)
        {
            if(Coliniar(p, q, this->points[i]) == 0)
            {
                r = this->points[i];
                r_id = i;
                break;
            }
        }
        //TODO: build the inner and outer faces, as well as other things in DCEL with p, q, r

        //TODO: allocate the points into buckets except p, q, r

    }

    ~DelaunayTriangulation()
    {
        for(int i = 0; i < this->points.size(); i ++)
        {
            delete(this->points[i]);
        }
        for(int i = 0; i < this->vertexs.size(); i ++)
        {
            delete(this->vertexs[i]);
        }
        for(int i = 0; i < this->half_edges.size(); i ++)
        {
            delete(this->half_edges[i]);
        }
        for(int i = 0; i < this->faces.size(); i ++)
        {
            delete(this->faces[i]);
        }
        this->points.clear();
        this->unassigned_points.clear();
        this->faces.clear();
        this->vertexs.clear();
        this->half_edges.clear();
    }
}









int main()
{
    vector<Point> points;
    points.clear();
    int n;
    scanf("%d", &n);
    for(int i = 1; i <= n; i ++)
    {
        long long x, y;
        scanf("%lld %lld", &x, &y);
        Point new_point = Point(i, x, y);
        points.push_back(new_point);
    }
    DelaunayTriangulation* dt = new DelaunayTriangulation(points);
    delete(dt);
    return 0;
}