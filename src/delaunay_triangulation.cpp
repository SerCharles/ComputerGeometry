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
    bool valid = 1; //whether the edge is valid
    HalfEdge* twin = NULL; //reference to the twin half edge
    Vertex* source = NULL; //reference to the source vertex
    Vertex* target = NULL;
    Face* incident_face = NULL; //reference to the left incident face
    HalfEdge* previous_halfedge = NULL; //reference to CCW previous half-edge
    HalfEdge* successor_halfedge = NULL; //reference to CCW next half-edge
    HalfEdge()
    {
        this->valid = 1;
    }
    ~HalfEdge(){}

};

class Face
{
public:
    bool valid = 1; //whether the face is valid
    HalfEdge* incident_halfedge = NULL; //the first half edge incident to the face from left
    vector<Point*> bucket_points; //the unassigned points of the bucket
    Face()
    {
        this->valid = 1;
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
    Point* p = pq->source->p;
    Point* q = qr->source->p;
    Point* r = rp->source->p;
    bool result = 0;
    if(ToLeft(p, q, point) && ToLeft(q, r, point) && ToLeft(r, p, point))
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
    d [Point*]: [one of the four points]

Returns:
    result [bool]: [inside: 1, not inside: 0]
*/
bool InCircle(Point* a, Point* b, Point* c, Point* d)
{
    double matrix[3][3];
    matrix[0][0] = a->x - d->x;
    matrix[0][1] = a->y - d->y;
    matrix[0][2] = (a->x - d->x) * (a->x - d->x) + (a->y - d->y) * (a->y - d->y);
    matrix[1][0] = b->x - d->x;
    matrix[1][1] = b->y - d->y;
    matrix[1][2] = (b->x - d->x) * (b->x - d->x) + (b->y - d->y) * (b->y - d->y);
    matrix[2][0] = c->x - d->x;
    matrix[2][1] = c->y - d->y;
    matrix[2][2] = (c->x - d->x) * (c->x - d->x) + (c->y - d->y) * (c->y - d->y);
    //threshold = 1e-10
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
    vector<Point*> points; //all the points
    vector<Point*> unassigned_points; //the unassigned points
    vector<Vertex*> all_vertexs; //all vertexs, used in recording and deletion
    vector<HalfEdge*> all_half_edges; //all half edges, used in recording and deletion
    vector<Face*> all_faces; //all the existed faces, used in recording and deletion
    int n_points; //the number of points, not including the for bounding box points

    DelaunayTriangulation(vector<Point>& points_)
    {
        //record all the points
        this->points.clear();
        this->unassigned_points.clear();
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
        this->all_vertexs.clear();
        this->all_half_edges.clear();
        this->all_faces.clear();
    }
    void InitTriangle(int size);
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
    pq->source = p;
    qr->source = q;
    rp->source = r;
    pq->target = q;
    qr->target = r;
    rp->target = p;
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
    pr->source = p;
    rs->source = r;
    sp->source = s;
    pr->target = r;
    rs->target = s;
    sp->target = p;
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
    HalfEdge* ab = abc->incident_halfedge;
    HalfEdge* bc = ab->successor_halfedge;
    HalfEdge* ca = bc->successor_halfedge;
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
    p->incident_halfedge = pa;
    pa->source = p;
    pb->source = p;
    pc->source = p;
    ap->source = a;
    bp->source = b;
    cp->source = c;
    pa->target = a;
    pb->target = b;
    pc->target = c;
    ap->target = p;
    bp->target = p;
    cp->target = p;
    pa->incident_face = pab;
    pb->incident_face = pbc;
    pc->incident_face = pca;
    ap->incident_face = pca;
    bp->incident_face = pab;
    cp->incident_face = pbc;
    ab->incident_face = pab;
    bc->incident_face = pbc;
    ca->incident_face = pca;
    pa->twin = ap;
    ap->twin = pa;
    pb->twin = bp;
    bp->twin = pb;
    pc->twin = cp;
    cp->twin = pc;
    pa->successor_halfedge = ab;
    ab->successor_halfedge = bp;
    bp->successor_halfedge = pa;
    pa->previous_halfedge = bp;
    ab->previous_halfedge = pa;
    bp->previous_halfedge = ab;
    pb->successor_halfedge = bc;
    bc->successor_halfedge = cp;
    cp->successor_halfedge = pb;
    pb->previous_halfedge = cp;
    bc->previous_halfedge = pb;
    cp->previous_halfedge = bc;
    pc->successor_halfedge = ca;
    ca->successor_halfedge = ap;
    ap->successor_halfedge = pc;
    pc->previous_halfedge = ap;
    ca->previous_halfedge = pc;
    ap->previous_halfedge = ca;
    pab->incident_halfedge = pa;
    pbc->incident_halfedge = pb;
    pca->incident_halfedge = pc;
    this->all_half_edges.push_back(pa);
    this->all_half_edges.push_back(pb);
    this->all_half_edges.push_back(pc);
    this->all_half_edges.push_back(ap);
    this->all_half_edges.push_back(bp);
    this->all_half_edges.push_back(cp);
    this->all_faces.push_back(pab);
    this->all_faces.push_back(pbc);
    this->all_faces.push_back(pca);


    //update bucket points
    for(int i = 0; i < abc->bucket_points.size(); i ++)
    {
        Point* the_point = abc->bucket_points[i];
        if(the_point->bucket == NULL)
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
    abc->valid = 0;
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
    HalfEdge* ax = ba->successor_halfedge;
    HalfEdge* xb = ax->successor_halfedge;
    HalfEdge* pa = ab->previous_halfedge;
    HalfEdge* bp = ab->successor_halfedge;
    Face* pab = pa->incident_face;
    Face* bax = ba->incident_face;
    HalfEdge* px = new HalfEdge();
    HalfEdge* xp = new HalfEdge();
    px->source = p;
    xp->source = x;
    px->target = x;
    xp->target = p;
    pa->successor_halfedge = ax;
    ax->successor_halfedge = xp;
    xp->successor_halfedge = pa;
    pa->previous_halfedge = xp;
    ax->previous_halfedge = pa;
    xp->previous_halfedge = ax;
    bp->successor_halfedge = px;
    px->successor_halfedge = xb;
    xb->successor_halfedge = bp;
    bp->previous_halfedge = xb;
    px->previous_halfedge = bp;
    xb->previous_halfedge = px;
    px->twin = xp;
    xp->twin = px;
    Face* pax = new Face();
    Face* bpx = new Face();
    pax->incident_halfedge = pa;
    bpx->incident_halfedge = bp;
    pa->incident_face = pax;
    ax->incident_face = pax;
    xp->incident_face = pax;
    bp->incident_face = bpx;
    px->incident_face = bpx;
    xb->incident_face = bpx;
    this->all_half_edges.push_back(px);
    this->all_half_edges.push_back(xp);
    this->all_faces.push_back(pax);
    this->all_faces.push_back(bpx); 

    //update buckets
    for(int i = 0; i < pab->bucket_points.size(); i ++)
    {
        if(InTriangle(pax, pab->bucket_points[i]))
        {
            pax->bucket_points.push_back(pab->bucket_points[i]);
            pab->bucket_points[i]->bucket = pax;
        }
        else 
        {
            bpx->bucket_points.push_back(pab->bucket_points[i]);
            pab->bucket_points[i]->bucket = bpx;
        }
    }
    for(int i = 0; i < bax->bucket_points.size(); i ++)
    {
        if(InTriangle(pax, bax->bucket_points[i]))
        {
            pax->bucket_points.push_back(bax->bucket_points[i]);
            bax->bucket_points[i]->bucket = pax;
        }
        else 
        {
            bpx->bucket_points.push_back(bax->bucket_points[i]);
            bax->bucket_points[i]->bucket = bpx;
        }
    }
    pab->bucket_points.clear();
    bax->bucket_points.clear();
    pab->valid = 0;
    bax->valid = 0;
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
    if(ab == NULL)
    {
        return NULL;
    }
    HalfEdge* ba = ab->twin;
    if(ba == NULL)
    {
        return NULL;
    }
    HalfEdge* xb = ba->previous_halfedge;
    if(xb == NULL)
    {
        return NULL;
    }
    Vertex* x = xb->source;
    return x;
}

/*
Iteratively insert one point to update the delaunay triangulation
*/
void DelaunayTriangulation::InsertOnePoint()
{
    //get the new point and the bucket faces
    Point* new_point = this->unassigned_points[unassigned_points.size() - 1];
    this->unassigned_points.pop_back();
    Face* original_face = new_point->bucket;
    new_point->bucket = NULL;
    HalfEdge* ab = original_face->incident_halfedge;
    HalfEdge* bc = ab->successor_halfedge;
    HalfEdge* ca = bc->successor_halfedge;
    Vertex* p = new Vertex(new_point);
    this->all_vertexs.push_back(p);
    Vertex* a = ab->source;
    Vertex* b = bc->source;
    Vertex* c = ca->source;

    //connect pa, pb, pc, form 3 new triangles
    this->ConnectTriangle(p, original_face);

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
        if(x == NULL)
        {
            continue;
        }
        if(InCircle(p->p, a->p, b->p, x->p))
        {
            HalfEdge* ba = ab->twin;
            HalfEdge* ax = ba->successor_halfedge;
            HalfEdge* xb = ax->successor_halfedge;
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
            if(s >= 0 && t >= 0)
            {
                result_edges.push_back(the_edge);
            }
        }
    }
    return result_edges;
}


struct UnionFind {
public:
    void clear(int n) {
        for (int i = 0; i <= n; i++) fa[i] = i;
    }
    void join(int x, int y) {
        int a = find(x), b = find(y);
        fa[a] = b;
    }
    int find(int x) { return (fa[x] == x) ? x : fa[x] = find(fa[x]); }
    int fa[100010];
};
UnionFind uf;
struct E2 {
    int x, y;
    double w;
    bool operator<(const E2& e) const { return w < e.w; }
};




int main()
{
    vector<Point> points;
    points.clear();
    int n;
    scanf("%d", &n);
    //shake the input to avoid coliniar, randomly shuffle the list
    mt19937 random_generator;
    random_generator.seed(14530529); //set seed to reproduce
    normal_distribution<double> nd(0, 1e-10); //a normal distribution noise
    for(int i = 1; i <= n; i ++)
    {
        int x, y;
        scanf("%d %d", &x, &y);
        double shaked_x = double(x) + nd(random_generator);
        double shaked_y = double(y) + nd(random_generator);
        Point new_point = Point(i, shaked_x, shaked_y);
        //Point new_point = Point(i, double(x), double(y));
        points.push_back(new_point);
    }
    shuffle(points.begin(), points.end(), random_generator);
    DelaunayTriangulation* dt = new DelaunayTriangulation(points);
    for(int i = 1; i <= n; i ++)
    {
        dt->InsertOnePoint();
    }
    vector<HalfEdge*> result_edges = dt->GetFinalEdges();

    
    std::vector<E2> edges;
    edges.clear();
    for(int i = 0; i < result_edges.size(); i ++)
    {
        Point* s = result_edges[i]->source->p;
        Point* t = result_edges[i]->target->p;
        int s_id = s->id;
        int t_id = t->id;
        if(s_id >= 0 && t_id >= 0)
        {
            double dist = sqrt((s->x - t->x) * (s->x - t->x) + (s->y - t->y) * (s->y - t->y));
            edges.push_back({s_id, t_id, dist});
        }
    }
    sort(edges.begin(), edges.end());
    double ans = 0;
    uf.clear(n);
    for (auto& e : edges) {
        int p = uf.find(e.x);
        int q = uf.find(e.y);
        if (p != q) {
            ans += e.w;
            uf.join(e.x, e.y);
        }
    }
    printf("%.6lf\n", ans);
    

    /*
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
    */
    //delete(dt);
    return 0;
}