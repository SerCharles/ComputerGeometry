#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <queue>
using namespace std;


//DCEL data structure
class Vertex;
class HalfEdge;
class Face;

class Point
{
public:
    long long id;
    double x;
    double y;
    bool used = 0; //used or not
    Face* bucket = NULL; //the bucket of the point
    Point(long long id_, double x_, double y_)
    {
        this->id = id_;
        this->x = x_;
        this->y = y_;
    }
    ~Point(){}
};

class Vertex
{
public:
    Point* p = NULL; //the coordinate of the vertex
    HalfEdge* half_edge = NULL; //reference to the first outgoing incident half-edge
    Vertex(Point* p_)
    {
        this->p = p_;
        this->half_edge = NULL;
    }
    ~Vertex(){}

};

class HalfEdge 
{
public:
    bool valid = 1; //whether the edge is valid
    HalfEdge* twin = NULL; //reference to the twin half edge
    Vertex* source = NULL; //reference to the source vertex
    Vertex* target = NULL;
    Face* face = NULL; //reference to the left incident face
    HalfEdge* next = NULL; //reference to CCW next half-edge
    HalfEdge()
    {
        this->valid = 1;
    }
    ~HalfEdge(){}

};

class Face
{
public:
    HalfEdge* half_edge = NULL; //the first half edge incident to the face from left
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
    if(inner_product >= 0)
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
    HalfEdge* pq = triangle->half_edge;
    HalfEdge* qr = pq->next;
    HalfEdge* rp = qr->next;
    Point* p = pq->source->p;
    Point* q = qr->source->p;
    Point* r = rp->source->p;
    bool result = 0;
    bool pq_result = ToLeft(p, q, point);
    bool qr_result = ToLeft(q, r, point);
    bool rp_result = ToLeft(r, p, point);
    if(pq_result == qr_result && qr_result == rp_result)
    {   
        result = 1;
    }
    return result;
}



/*
Calculate the determinant of a 3 * 3 matrix

Args:
    matrix [double array] [3 * 3]: [the input matrix]

Returns:
    result [double]: [the determinant of the matrix]
*/
double Determinant(double matrix[3][3])
{
    double result = 0;
    result += matrix[0][0] * matrix[1][1] * matrix[2][2];
    result += matrix[0][1] * matrix[1][2] * matrix[2][0];
    result += matrix[0][2] * matrix[1][0] * matrix[2][1];
    result -= matrix[2][0] * matrix[1][1] * matrix[0][2];
    result -= matrix[2][1] * matrix[1][2] * matrix[0][0];
    result -= matrix[2][2] * matrix[1][0] * matrix[0][1];
    return result;
}

/*
Judge whether a point lies in the circumcircle of the three other points

Args:
    a [Point*]: [one of the four points]
    b [Point*]: [one of the four points]
    c [Point*]: [one of the four points]
    p [Point*]: [one of the four points]

Returns:
    result [bool]: [inside: 1, not inside: 0]
*/

bool InCircle(Point* a, Point* b, Point* c, Point* p)
{
    double matrix[3][3];
    matrix[0][0] = a->x - p->x;
    matrix[0][1] = a->y - p->y;
    matrix[0][2] = (a->x - p->x) * (a->x - p->x) + (a->y - p->y) * (a->y - p->y);
    //matrix[0][2] = a->x * a->x + a->y * a->y - p->x * p->x - p->y * p->y;
    matrix[1][0] = b->x - p->x;
    matrix[1][1] = b->y - p->y;
    //matrix[1][2] = b->x * b->x + b->y * b->y - p->x * p->x - p->y * p->y;
    matrix[1][2] = (b->x - p->x) * (b->x - p->x) + (b->y - p->y) * (b->y - p->y);
    matrix[2][0] = c->x - p->x;
    matrix[2][1] = c->y - p->y;
    //matrix[1][2] = c->x * c->x + c->y * c->y - p->x * p->x - p->y * p->y;
    matrix[2][2] = (c->x - p->x) * (c->x - p->x) + (c->y - p->y) * (c->y - p->y);
    double determinant = Determinant(matrix);
    bool result = 0;
    if(determinant >= 0)
    {
        result = 1;
    }
    return result;
}


//the main class
class DelaunayTriangulation
{
public:
    vector<Point*> points; //all the unassigned points
    vector<HalfEdge*> all_half_edges; //all half edges, used in recording and deletion

    DelaunayTriangulation(vector<Point>& points_)
    {
        //record all the points
        this->points.clear();
        this->all_half_edges.clear();
        for(int i = 0; i < points_.size(); i ++)
        {
            Point* new_point = new Point(points_[i].id, points_[i].x, points_[i].y);
            this->points.push_back(new_point);
        }

        //set the first two triangles to be very big, covering all possible points
        this->InitTriangle(1e8);
    }
    void InitTriangle(double size);
    void InsertOnePoint();
    void ConnectTriangle(Vertex* p, Face* original_face);
    void FlipEdge(Vertex* p, HalfEdge* ab, Vertex* x);
    Vertex* RightSite(HalfEdge* ab);
    vector<HalfEdge*> GetFinalEdges();
};

/*
Build the initial two triangles, allocate the vertexs, half edges and faces, and buckets

Args:
    size [int]: [the range of the scene]

*/
void DelaunayTriangulation::InitTriangle(double size)
{
    //init the points、vertexs and faces
    Point* point_p = new Point(-1, 3 * size, size);
    Point* point_q = new Point(-1, -3 * size, size);
    Point* point_r = new Point(-1, 0, -2 * size);
    Vertex* p = new Vertex(point_p);
    Vertex* q = new Vertex(point_q);
    Vertex* r = new Vertex(point_r);
    Face* pqr = new Face();

    //build the half edges
    HalfEdge* pq = new HalfEdge();
    HalfEdge* qr = new HalfEdge();
    HalfEdge* rp = new HalfEdge();
    pqr->half_edge = pq;
    pq->face = pqr;
    qr->face = pqr;
    rp->face = pqr;
    pq->source = p;
    pq->target = q;
    qr->source = q;
    qr->target = r;
    rp->source = r;
    rp->target = p;
    pq->next = qr;
    qr->next = rp;
    rp->next = pq;
    
    
    //allocate the unassigned points into buckets
    for(int i = 0; i < this->points.size(); i ++)
    {
        Point* point = this->points[i];
        point->bucket = pqr;
        pqr->bucket_points.push_back(point);
    }

    //put the vertexs、
    point_p->used = 1;
    point_q->used = 1;
    point_r->used = 1;
    this->all_half_edges.push_back(pq);
    this->all_half_edges.push_back(qr);
    this->all_half_edges.push_back(rp);
}

/*
When adding a new vertex p to a triangle abc, update the DCEL data structure

Args:
    p [Vertex*]: [the vertex to be inserted to the triangle abc]
    original_face [Face*]: [the original triangle abc]
*/
void DelaunayTriangulation::ConnectTriangle(Vertex* p, Face* original_face)
{
    Face* abc = original_face;
    HalfEdge* ab = abc->half_edge;
    HalfEdge* bc = ab->next;
    HalfEdge* ca = bc->next;
    Vertex* a = ab->source;
    Vertex* b = bc->source;
    Vertex* c = ca->source;
    HalfEdge* pa = new HalfEdge();
    HalfEdge* pb = new HalfEdge();
    HalfEdge* pc = new HalfEdge();
    HalfEdge* ap = new HalfEdge();
    HalfEdge* bp = new HalfEdge();
    HalfEdge* cp = new HalfEdge();
    Face* pab = new Face();
    Face* pbc = new Face();
    Face* pca = new Face();
    p->half_edge = pa;
    pa->source = p;
    pa->target = a;
    pb->source = p;
    pb->target = b;
    pc->source = p;
    pc->target = c;

    ap->source = a;
    ap->target = p;
    bp->source = b;
    bp->target = p;
    cp->source = c;
    cp->target = p;

    pa->face = pab;
    pb->face = pbc;
    pc->face = pca;
    ap->face = pca;
    bp->face = pab;
    cp->face = pbc;
    ab->face = pab;
    bc->face = pbc;
    ca->face = pca;

    pa->twin = ap;
    ap->twin = pa;
    pb->twin = bp;
    bp->twin = pb;
    pc->twin = cp;
    cp->twin = pc;

    pa->next = ab;
    ab->next = bp;
    bp->next = pa;

    pb->next = bc;
    bc->next = cp;
    cp->next = pb;

    pc->next = ca;
    ca->next = ap;
    ap->next = pc;

    pab->half_edge = pa;
    pbc->half_edge = pb;
    pca->half_edge = pc;

    this->all_half_edges.push_back(pa);
    this->all_half_edges.push_back(pb);
    this->all_half_edges.push_back(pc);
    this->all_half_edges.push_back(ap);
    this->all_half_edges.push_back(bp);
    this->all_half_edges.push_back(cp);

    //update bucket points
    for(int i = 0; i < abc->bucket_points.size(); i ++)
    {
        Point* the_point = abc->bucket_points[i];
        if(the_point->used == 1)
        {
            continue;
        }
        if(InTriangle(pab, the_point))
        {
            pab->bucket_points.push_back(the_point);
            the_point->bucket = pab;
        }
        else if(InTriangle(pbc, the_point))
        {
            pbc->bucket_points.push_back(the_point);
            the_point->bucket = pbc;
        }
        else 
        {
            pca->bucket_points.push_back(the_point);
            the_point->bucket = pca;
        }
    }
    abc->bucket_points.clear();
    delete(abc);
}

/*
Flip Edge:Remove the edge ab and ba, add edge px and xp, renew the triangles

Args:
    p [Vertex*]: [Vertex p]
    ab [HalfEdge*]: [Half Edge ab]
    x [Vertex*]: [Vertex x]
*/
void DelaunayTriangulation::FlipEdge(Vertex* p, HalfEdge* ab, Vertex* x)
{
    HalfEdge* ba = ab->twin;
    HalfEdge* ax = ba->next;
    HalfEdge* xb = ax->next;
    HalfEdge* bp = ab->next;
    HalfEdge* pa = bp->next;
    Face* pab = pa->face;
    Face* bax = ba->face;
    HalfEdge* px = new HalfEdge();
    HalfEdge* xp = new HalfEdge();

    px->source = p;
    px->target = x;
    xp->source = x;
    xp->target = p;

    pa->next = ax;
    ax->next = xp;
    xp->next = pa;

    bp->next = px;
    px->next = xb;
    xb->next = bp;

    px->twin = xp;
    xp->twin = px;
    Face* pax = new Face();
    Face* bpx = new Face();
    pax->half_edge = xp;
    bpx->half_edge = px;
    pa->face = pax;
    ax->face = pax;
    xp->face = pax;
    bp->face = bpx;
    px->face = bpx;
    xb->face = bpx;
    this->all_half_edges.push_back(px);
    this->all_half_edges.push_back(xp);

    //update buckets
    for(int i = 0; i < pab->bucket_points.size(); i ++)
    {
        Point* the_point = pab->bucket_points[i];
        if(the_point->used == 1)
        {
            continue;
        }
        if(InTriangle(pax, the_point))
        {
            pax->bucket_points.push_back(the_point);
            the_point->bucket = pax;
        }
        else 
        {
            bpx->bucket_points.push_back(the_point);
            the_point->bucket = bpx;
        }
    }
    for(int i = 0; i < bax->bucket_points.size(); i ++)
    {
        Point* the_point = bax->bucket_points[i];
        if(the_point->used == 1)
        {
            continue;
        }
        if(InTriangle(pax, the_point))
        {
            pax->bucket_points.push_back(the_point);
            the_point->bucket = pax;
        }
        else 
        {
            bpx->bucket_points.push_back(the_point);
            the_point->bucket = bpx;
        }
    }
    pab->bucket_points.clear();
    bax->bucket_points.clear();
    delete(pab);
    delete(bax);
    ab->valid = 0;
    ba->valid = 0;
}

/*
Get the right site of a half edge ab

Args:
    ab [HalfEdge*]: [the half edge to be searched]

Returns:
    x [Vertex*]: [the right site of ab]
*/
Vertex* DelaunayTriangulation::RightSite(HalfEdge* ab)
{
    HalfEdge* ba = ab->twin;
    if(ba == NULL || ba->valid == 0)
    {
        return NULL;
    }
    HalfEdge* ax = ba->next;
    
    if(ax == NULL || ax->valid == 0)
    {
        return NULL;
    }
    Vertex* x = ax->target;
    return x;
}

/*
Iteratively insert one point to update the delaunay triangulation
*/
void DelaunayTriangulation::InsertOnePoint()
{
    //get the new point and the bucket faces
    Point* new_point = this->points[this->points.size() - 1];
    this->points.pop_back();
    Face* abc = new_point->bucket;
    new_point->used = 1;
    Vertex* p = new Vertex(new_point);
    HalfEdge* ab = abc->half_edge;
    HalfEdge* bc = ab->next;
    HalfEdge* ca = bc->next;
    Vertex* a = ab->source;
    Vertex* b = bc->source;
    Vertex* c = ca->source;

    //connect pa, pb, pc, form 3 new triangles
    this->ConnectTriangle(p, abc);


    //Swap Test
    queue<HalfEdge*> Q;
    while(!Q.empty())
    {
        Q.pop();
    }
    Q.push(ab);
    Q.push(bc);
    Q.push(ca);
    while(!Q.empty())
    {
        HalfEdge* ab = Q.front();
        Vertex* a = ab->source;
        Vertex* b = ab->target;
        Q.pop();
        Vertex* x = this->RightSite(ab);
        if(x == NULL || x == p)
        {
            continue;
        }
        if(InCircle(p->p, a->p, b->p, x->p))
        {
            HalfEdge* ba = ab->twin;
            HalfEdge* ax = ba->next;
            HalfEdge* xb = ax->next;
            Q.push(ax);
            Q.push(xb);
            this->FlipEdge(p, ab, x);
        }
    }
}

/*
Get all the final valid edges of the result, only leaving one of the twins, remove those of the bounding box

Returns:
    result_edges [vector<HalfEdge*>]: [the result edges]
*/
vector<HalfEdge*> DelaunayTriangulation::GetFinalEdges()
{
    vector<HalfEdge*> result_edges;
    result_edges.clear();
    for(int i = 0; i < this->all_half_edges.size(); i ++)
    {
        HalfEdge* the_edge = this->all_half_edges[i];
        if(the_edge->valid == 1)
        {
            if(the_edge->twin != NULL)
            {
                the_edge->twin->valid = 0;
            }
            long long s = the_edge->source->p->id;
            long long t = the_edge->target->p->id;
            if(s >= 0 && t >= 0 && s != t)
            {
                result_edges.push_back(the_edge);
            }
        }
    }
    return result_edges;
}


int main()
{
    vector<Point> points;
    points.clear();
    int n;
    scanf("%d", &n);
    for(int i = 1; i <= n; i ++)
    {
        int x, y;
        scanf("%d %d", &x, &y);
        Point new_point = Point(i, double(x), double(y));
        points.push_back(new_point);
    }

    //randomly shuffle the list
    mt19937 random_generator;
    random_generator.seed(14530529); //set seed to reproduce
    shuffle(points.begin(), points.end(), random_generator);

    DelaunayTriangulation* dt = new DelaunayTriangulation(points);
    for(int i = 1; i <= n; i ++)
    {
        dt->InsertOnePoint();
    }

    vector<HalfEdge*> result_edges = dt->GetFinalEdges();
    long long result = 0;
    for(int i = 0; i < result_edges.size(); i ++)
    {
        long long s = result_edges[i]->source->p->id;
        long long t = result_edges[i]->target->p->id;
        result += s;
        result += t;
    }
    result = result % (result_edges.size() + 1);
    printf("%lld", result);
    return 0;
}