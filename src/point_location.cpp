#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
#define POINT 1
#define LINE 2
#define TRAP 3

/*
Definition of a point
Args:
    x [double]: [the x coordinate of the point]
    y [double]: [the y coordinate of the point]
*/
class Point
{
public:
    double x;
    double y;
    Point() 
    {
        this->x = 0;
        this->y = 0;
    }
    Point(double a, double b)
    {
        this->x = a;
        this->y = b;
    }
    ~Point()
    {

    }

};

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
Definition of a line
Args:
    id [long long]: [the id of the line]
    s [Point]: [the starting point of the line, which is the left or down end]
    t [Point]: [the ending point of the line, which is the right or up end]
*/
class Line
{
public:
    long long id;
    Point* s;
    Point* t;
    bool vertical; //whether the line is vertical or not
    Line()
    {
        this->id = 0;
        this->vertical = 1;
    }
    Line(long long id_, Point* s_, Point* t_)
    {
        this->id = id_;
        if(s_->x != t_->x)
        {
            this->vertical = 0;
            if(s_->x < t_->x)
            {
                this->s = s_;
                this->t = t_;
            }
            else 
            {
                this->s = t_;
                this->t = s_;
            }
        }
        else 
        {
            this->vertical = 1;
            if(s_->y < t_->y)
            {
                this->s = s_;
                this->t = t_;
            }
            else 
            {
                this->s = t_;
                this->t = s_;
            }
        }
    }
    ~Line()
    {
        delete(this->s);
        delete(this->t);
    }
};

//Used in sorting for special judge
bool CompareSX(Line* a, Line* b)
{
    return a->s->x < b->s->x;
}
bool CompareTX(Line* a, Line* b)
{
    return a->t->x < b->t->x;
}


/*
Definition of a trapezoid
Args:
    left [double]: [the left x of the trapezoid]
    right [double]: [the right x of the trapezoid]
    up [Line]: [the up line of the trapezoid]
    down [Line]: [the down line of the trapezoid]
*/
class Trap
{
public:
    double left;
    double right;
    Line* up = NULL;
    Line* down = NULL;
    Trap()
    {
        this->left = 0;
        this->right = 0;
        this->up = NULL;
        this->down = NULL;
    }
    Trap(double left_, double right_, Line* up_, Line* down_)
    {
        this->left = left_;
        this->right = right_;
        this->up = up_;
        this->down = down_;
    }
    ~Trap()
    {

    }
};



/*
The father class of all Trapezoid nodes
*/
class Node
{
public:
    int type; //the type of each node, 1:point 2:line 3:trap
    Node* left_son = NULL;
    Node* right_son = NULL;
    vector<Node*> fathers;
    Node()
    {
        this->type = 0;
        this->left_son = NULL; 
        this->right_son = NULL; 
        this->fathers.clear();
    }
    ~Node()
    {
        this->fathers.clear();
        delete(this->left_son);
        delete(this->right_son);
    }
    

};

/*
The point node
*/
class PointNode : public Node 
{
public:
    Point* point = NULL; //the point of the node
    PointNode()
    {
        this->type = POINT;
        this->point = NULL;
    }
    PointNode(Point* p)
    {
        this->type = POINT;
        this->point = p;
    }
    ~PointNode()
    {

    }
};

/*
The line node
*/
class LineNode : public Node 
{
public:
    Line* line = NULL;
    LineNode()
    {
        this->type = LINE;
        this->line = NULL;
    }
    LineNode(Line* l)
    {
        this->type = LINE;
        this->line = l;
    }
    ~LineNode()
    {

    }
};

/*
The trap node
*/
class TrapNode : public Node 
{
public:
    Trap* trap = NULL;
    TrapNode* left_up = NULL; //the left upper trap of the node
    TrapNode* left_down = NULL; //the left down trap of the node
    TrapNode* right_up = NULL; //the right upper trap of the node
    TrapNode* right_down = NULL; //the right down trap of the node
    TrapNode()
    {
        this->type = TRAP;
        this->trap = NULL;
        this->left_up = NULL;
        this->left_down = NULL;
        this->right_up = NULL;
        this->right_down = NULL;
    }
    TrapNode(Trap* t)
    {
        this->type = TRAP;
        this->trap = t;
        this->left_up = NULL;
        this->left_down = NULL;
        this->right_up = NULL;
        this->right_down = NULL;
    }
    ~TrapNode()
    {
        delete(trap);
    }
};

//the definition of the search structure
class TrapezoidalTree
{
public:
    //the range of the data, used in initialization
    double min_x = -1000000;
    double max_x = 1000000;
    double min_y = -1000000;
    double max_y = 1000000; 
    Line* up_line = NULL;
    Line* down_line = NULL;

    vector<Line*> lines; // all the lines for the trapezoidal map

    Node* root = NULL; //the root of the tree

    TrapezoidalTree(vector<Line*>& lines_)
    {
        //record the lines
        for(int i = 0; i < lines_.size(); i ++)
        {
            this->lines.push_back(lines_[i]);
        }

        //init the framework: the root node is a trap node of the bounding box of the tree
        Point* up_left = new Point(min_x, max_y);
        Point* up_right = new Point(max_x, max_y);
        Point* down_left = new Point(min_x, min_y);
        Point* down_right = new Point(max_x, min_y);
        this->up_line = new Line(0, up_left, up_right);
        this->down_line = new Line(0, down_left, down_right);
        Trap* main_trap = new Trap(min_x, max_x, this->up_line, this->down_line);
        this->root = new TrapNode(main_trap);
    } 

    ~TrapezoidalTree()
    {
        delete(this->root);
        delete(this->up_line);
        delete(this->down_line);
    }

    void BuildTree();
    TrapNode* SearchPoint(Point* p);
    void InsertLine(TrapNode* trap, Line* new_line);

};

/*
Search the place of a point in the subtree of the trapezoidal tree

Args:
    node [Node*]: [the root of the subtree]
    p [Point*]: [the point to be searched]

Returns:
    place [TrapNode*]: [the trap node to be searched]
*/
TrapNode* SearchPoint(Node* node, Point* p)
{
    
    if(node->type == TRAP)
    {
        //trap node, leaf node
        return (TrapNode*)node;
    }
    else if(node->type == POINT)
    {
        //point node
        //< left, >= right: buggy, need to special judge the right point of each lines by binary search
        double point_x = (PointNode*)node->point->x;
        if(p->x < point_x)
        {
            return this->SearchPoint(node->left_son, p);
        }
        else 
        {
            return this->SearchPoint(node->right_son, p);
        }
    }
    else 
    {
        //line node
        Line* line = (LineNode*)node->line;
        //strictly up: left, others: right
        double area2 = Area2(line->s, line->t, p);
        if(area2 > 0)
        {
            return this->SearchPoint(node->left_son, p);
        }
        else 
        {
            return this->SearchPoint(node->right_son, p);
        }
    }
}

/*
Build the trapezoidal tree by the RIC algorithm
*/
void TrapezoidalTree::BuildTree()
{
    for(int i = 0; i < lines.size(); i ++)
    {
        Line* l = lines[i];
        Point* p = l->s;
        TrapNode* target_trap = this->SearchPoint(p);
        this->InsertLine(target_trap, l);
    }
}

/*
Insert a new line into the trapezoidal tree

Args:
    trap [TrapNode*]: [the trap containing the line's left point]
    new_line [Line*]: [the new line to be inserted]
*/
void TrapezoidalTree::InsertLine(TrapNode* trap, Line* new_line)
{
    //find all relavant trapezoidal nodes



}

/*
Special judge the lines, when the left or right point lies straight ahead of the query point, solving the bugs of border and vertical lines

Args:
    query_point [Point*]: [the point to be queried]
    lines_sx [vector<Line*>]: [the line list sorted by the x or the starting point(left or down)]
    lines_tx [vector<Line*>]: [the line list sorted by the x or the ending point(right or up)]

Returns:
    best_line [Line*]: [the best line of the special search, NULL if nothing]
    best_dist [double]: [the min distance of the special search, >= 0]
*/
Line* SpecialJudge(Point* query_point, vector<Line*>& lines_sx, vector<Line*>& lines_tx, double* best_dist)
{
    Line* best_line = NULL;
    best_dist = 2e6;
    double target_x = query_point->x;
    double target_y = query_point->y;

    //search in sx
    int start = 0;
    int end = lines_sx.size() - 1;
    while(start <= end)
    {
        int middle = (start + end) / 2;
        double current_x = lines_sx[middle]->s->x;
        if(current_x < target_x)
        {
            start = middle + 1;
        }
        else 
        {
            end = middle - 1;
        }
    }
    int current_place = start;
    while(lines_sx[current_place]->s->x == target_x)
    {
        double current_y = lines_sx[current_place]->s->y;
        double dist_y = current_y - target_y;

        //judge vertical lines
        if(lines_sx[current_place]->vertical == 1)
        {
            if(current_y >= lines_sx[current_place]->s->y && current_y <= lines_sx[current_place]->t->y)
            {
                dist_y = 0;
            }
        }

        if(dist_y >= 0 && dist_y < best_dist)
        {
            best_dist = dist_y;
            best_line = lines_sx[current_place];
        }
        current_place += 1;
    }

    //search in tx
    start = 0;
    end = lines_tx.size() - 1;
    while(start <= end)
    {
        int middle = (start + end) / 2;
        double current_x = lines_tx[middle]->t->x;
        if(current_x < target_x)
        {
            start = middle + 1;
        }
        else 
        {
            end = middle - 1;
        }
    }
    int current_place = start;
    while(lines_tx[current_place]->t->x == target_x)
    {
        double current_y = lines_tx[current_place]->t->y;
        double dist_y = current_y - target_y;

        //judge vertical lines
        if(lines_sx[current_place]->vertical == 1)
        {
            if(current_y >= lines_sx[current_place]->s->y && current_y <= lines_sx[current_place]->t->y)
            {
                dist_y = 0;
            }
        }

        if(dist_y >= 0 && dist_y < best_dist)
        {
            best_dist = dist_y;
            best_line = lines_tx[current_place];
        }
        current_place += 1;
    }
    return best_line;
}






int main()
{
    //input
    int n, m;
    scanf("%d %d", &n, &m);
    vector<Line*> normal_lines; //non vertical lines
    vector<Line*> all_lines_sx; //all lines sort by s.x, used in special judge
    vector<Line*> all_lines_tx; //all lines sort by t.x, used in special judge
    vector<Point*> query_points;
    normal_lines.clear();
    all_lines.clear();
    query_points.clear();
    for(int i = 1; i <= n; i ++)
    {
        int sx, sy, tx, ty;
        scanf("%d %d %d %d", &sx, &sy, &tx, &ty);
        Point* s = new Point(double(sx), double(sy));
        Point* t = new Point(double(tx), double(ty));
        Line* new_line = new Line(i, s, t);
        if(new_line->vertical == 0)
        {
            normal_lines.push_back(new_line);
        }
        all_lines_sx.push_back(new_line);
        all_lines_tx.push_back(new_line);
    }
    for(int i = 0; i < m; i ++)
    {
        int x, y;
        scanf("%d %d", &x, &y);
        Point* s = new Point(double(x), double(y));
        query_points.push_back(s);
    }
    sort(all_lines_sx.begin(), all_lines_sx.end(), CompareSX);
    sort(all_lines_tx.begin(), all_lines_tx.end(), CompareTX);


    //main
    TrapezoidalTree* the_tree = new TrapezoidalTree(normal_lines);

    //delete
    delete(the_tree);
    for(int i = 0; i < normal_lines.size(); i ++)
    {
        delete(normal_lines[i]);
    }
    for(int i = 0; i < all_lines_sx.size(); i ++)
    {
        delete(all_lines_sx[i]);
    }
    for(int i = 0; i < all_lines_tx.size(); i ++)
    {
        delete(all_lines_tx[i]);
    }
    for(int i = 0; i < query_points.size(); i ++)
    {
        delete(query_points[i]);
    }
    normal_lines.clear();
    all_lines.clear();
    query_points.clear();
    return 0;
}