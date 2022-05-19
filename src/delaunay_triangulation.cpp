#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;


//DCEL data structure
class Vertex;
class HalfEdge;
class Face;

class Point
{
public:
    int id;
    double x;
    double y;
    Face* bucket = NULL; //the bucket of the point
    Point(int id_, double x_, double y_)
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
};

//useful functions
/*
Get the signed area of a triangle defined by three points

Args:
    p [Point*]: [one of the three points of the triangle]
    q [Point*]: [one of the three points of the triangle]
    s [Point*]: [one of the three points of the triangle]
    
Returns:
    area [double]: [the signed area * 2 of the triangle, it is positive if s is at the left of the line pq]
*/
double Area2(Point* p, Point* q, Point* s)
{
    double area = p->x * q->y - p->y * q->x + q->x * s->y - q->y * s->x + s->x * p->y - s->y * p->x;
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
    double inner_product = (p->x - s->x) * (s->x - q->x) + (p->y - s->y) * (s->y - q->y);
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

/*
The In triangle test, testing whether a point is inside a triangle 

Args:
    triangle [Face*]: [the triangle]
    point [Point*]: [the point to be tested]
    
Returns:
    result [bool]: [1 if the point is inside the triangle]
*/
bool InTriangle(Face* triangle, Point* point)
{
    HalfEdge* pq = triangle->incident_halfedge;
    HalfEdge* qr = pq->successor_halfedge;
    HalfEdge* rp = qr->successor_halfedge;
    Point* p = pq->origin->p;
    Point* q = qr->origin->p;
    Point* r = rp->origin->p;
    bool result = 0;
    if(ToLeft(p, q, point) && ToLeft(q, r, point) && ToLeft(r, p, point))
    {   
        result = 1;
    }
    return result;
}



//the main class
class DelaunayTriangulation
{
public:
    vector<Point*> points; //all the points
    vector<Point*> unassigned_points; //the unassigned points
    vector<Face*> faces; //the used faces
    vector<Vertex*> all_vertexs; //all vertexs, used in recording and deletion
    vector<HalfEdge*> all_half_edges; //all half edges, used in recording and deletion
    vector<Face*> all_faces; //all the existed faces, used in recording and deletion
    int n_points; //the number of points, not including the for bounding box points

    DelaunayTriangulation(vector<Point>& points_)
    {
        //record all the points
        this->points.clear();
        this->unassigned_points.clear();
        this->faces.clear();
        this->all_vertexs.clear();
        this->all_half_edges.clear();
        this->all_faces.clear();
        this->n_points = points_.size();
        for(int i = 0; i < n_points; i ++)
        {
            Point* new_point = new Point(points_[i].id, points_[i].x, points_[i].y);
            this->points.push_back(new_point);
        }

        //set the first two triangles to be very big, covering all possible points
        this->InitTriangle(1e7);
    }

    ~DelaunayTriangulation()
    {
        for(int i = 0; i < this->points.size(); i ++)
        {
            delete(this->points[i]);
        }
        for(int i = 0; i < this->all_vertexs.size(); i ++)
        {
            delete(this->all_vertexs[i]);
        }
        for(int i = 0; i < this->all_half_edges.size(); i ++)
        {
            delete(this->all_half_edges[i]);
        }
        for(int i = 0; i < this->all_faces.size(); i ++)
        {
            delete(this->all_faces[i]);
        }
        this->points.clear();
        this->unassigned_points.clear();
        this->faces.clear();
        this->all_vertexs.clear();
        this->all_half_edges.clear();
        this->all_faces.clear();
    }
    void InitTriangle(int size);
};

/*
Build the initial two triangles, allocate the vertexs, half edges and faces, and buckets

Args:
    size [int]: [the range of the scene]

*/
void DelaunayTriangulation::InitTriangle(int size)
{
    //init the points、vertexs and faces
    Point* point_p = new Point(-1, -1e7, -1e7);
    Point* point_q = new Point(-1, 1e7, -1e7);
    Point* point_r = new Point(-1, 1e7, 1e7);
    Point* point_s = new Point(-1, -1e7, 1e7);
    Vertex* p = new Vertex(point_p);
    Vertex* q = new Vertex(point_q);
    Vertex* r = new Vertex(point_r);
    Vertex* s = new Vertex(point_s);
    Face* face_pqr = new Face();
    Face* face_prs = new Face();

    //build the half edges
    HalfEdge* pq = new HalfEdge();
    HalfEdge* qr = new HalfEdge();
    HalfEdge* rp = new HalfEdge();
    face_pqr->incident_halfedge = pq;
    pq->incident_face = face_pqr;
    qr->incident_face = face_pqr;
    rp->incident_face = face_pqr;
    pq->origin = p;
    qr->origin = q;
    rp->origin = r;
    pq->previous_halfedge = rp;
    qr->previous_halfedge = pq;
    rp->previous_halfedge = qr;
    pq->successor_halfedge = qr;
    qr->successor_halfedge = rp;
    rp->successor_halfedge = pq;
    HalfEdge* pr = new HalfEdge();
    HalfEdge* rs = new HalfEdge();
    HalfEdge* sp = new HalfEdge();
    face_prs->incident_halfedge = pr;
    pr->incident_face = face_prs;
    rs->incident_face = face_prs;
    sp->incident_face = face_prs;
    pr->origin = p;
    rs->origin = r;
    sp->origin = s;
    pr->previous_halfedge = sp;
    rs->previous_halfedge = pr;
    sp->previous_halfedge = rs;
    pr->successor_halfedge = rs;
    rs->successor_halfedge = sp;
    sp->successor_halfedge = pr;
    pr->twin = rp;
    rp->twin = pr;
    


    //allocate the unassigned points into buckets
    for(int i = 0; i < this->points.size(); i ++)
    {
        Point* point = this->points[i];
        if(InTriangle(face_pqr, point))
        {
            point->bucket = face_pqr;
            face_pqr->bucket_points.push_back(point);
        }
        else 
        {
            point->bucket = face_prs;
            face_prs->bucket_points.push_back(point);
        }
        this->unassigned_points.push_back(point);
    }

    //put the vertexs、
    this->points.push_back(point_p);
    this->points.push_back(point_q);
    this->points.push_back(point_r);
    this->points.push_back(point_s);
    this->all_vertexs.push_back(p);
    this->all_vertexs.push_back(q);
    this->all_vertexs.push_back(r);
    this->all_vertexs.push_back(s);
    this->all_half_edges.push_back(pq);
    this->all_half_edges.push_back(qr);
    this->all_half_edges.push_back(rp);
    this->all_half_edges.push_back(pr);
    this->all_half_edges.push_back(rs);
    this->all_half_edges.push_back(sp);
    this->all_faces.push_back(face_pqr);
    this->all_faces.push_back(face_prs);
    this->faces.push_back(face_pqr);
    this->faces.push_back(face_prs);
}








int main()
{
    vector<Point> points;
    points.clear();
    int n;
    scanf("%d", &n);
    //shake the input to avoid coliniar
    mt19937 rng(random_device{}());
    normal_distribution<double> nd(0, 1e-10);
    for(int i = 1; i <= n; i ++)
    {
        int x, y;
        scanf("%d %d", &x, &y);
        double shaked_x = double(x) + nd(rng);
        double shaked_y = double(y) + nd(rng);
        Point new_point = Point(i, shaked_x, shaked_y);
        points.push_back(new_point);
    }
    DelaunayTriangulation* dt = new DelaunayTriangulation(points);
    delete(dt);
    return 0;
}