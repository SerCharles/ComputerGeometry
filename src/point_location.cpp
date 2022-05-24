#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <list>
#include <algorithm>
#include <random>
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

    /*
    Get the intersection y of the line and a vertical line x = t
    
    Args:
        x [double]: [the line x = t]
    
    Returns:
        y [double]: [the y of the intersection point]
    */
    double GetIntersectionY(double x)
    {
        double y = this->s->y + (x - this->s->x) / (this->t->x - this->s->x) * (this->t->y - this->s->y);
        return y;
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
        if(this->left_son != NULL && this->left_son->type != TRAP)
        {
            delete(this->left_son);
        }
        if(this->right_son != NULL && this->right_son->type != TRAP)
        {
            delete(this->right_son);
        }
        this->fathers.clear();
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
        if(this->trap != NULL)
        {
            delete(this->trap);
        }
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
        this->down_line = new Line(-1, down_left, down_right);
        Trap* main_trap = new Trap(min_x, max_x, this->up_line, this->down_line);
        this->root = new TrapNode(main_trap);
    } 

    void BuildTree();
    TrapNode* SearchPoint(Node* node, Point* p);
    void InsertLine(TrapNode* trap, Line* new_line);
    void UpdateFatherSon(Node* old_node, Node* new_node);
    void UpdateNeighborTopology(TrapNode* old_node, TrapNode* new_node, TrapNode* neighbor_node);
    bool JudgeTrapNodeMerge(TrapNode* left, TrapNode* right);
    void MergeTrapNode(TrapNode* left, TrapNode* right);
};

/*
Search the place of a point in the subtree of the trapezoidal tree

Args:
    node [Node*]: [the root of the subtree]
    p [Point*]: [the point to be searched]

Returns:
    place [TrapNode*]: [the trap node to be searched]
*/
TrapNode* TrapezoidalTree::SearchPoint(Node* node, Point* p)
{
    
    if(node->type == TRAP)
    {
        //trap node, leaf node
        return (TrapNode*)node;
    }
    else if(node->type == POINT)
    {
        //point node
        //< left, >= right: buggy for border points, need to special judge by binary search
        double point_x = ((PointNode*)node)->point->x;
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
        Line* line = ((LineNode*)node)->line;
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
        TrapNode* target_trap = this->SearchPoint(this->root, p);
        this->InsertLine(target_trap, l);
    }
}

/*
Update the topology of the neighbor of the inserted node

Args:
    old_node [TrapNode*]: [the old node]
    new_node [TrapNode*]: [the inserted new node]
    neighbor_node [TrapNode*]: [the neighbor node of the old node to be updated]
*/
void TrapezoidalTree::UpdateNeighborTopology(TrapNode* old_node, TrapNode* new_node, TrapNode* neighbor_node)
{
    if (neighbor_node == NULL) 
    {
        return;
    }
    if (neighbor_node->left_up == old_node)
    {
        neighbor_node->left_up = new_node;
    }
    if (neighbor_node->left_down == old_node)
    {
        neighbor_node->left_down = new_node;
    }
    if (neighbor_node->right_up == old_node)
    {
        neighbor_node->right_up = new_node;
    }
    if (neighbor_node->right_down == old_node)
    {
        neighbor_node->right_down = new_node;
    }
}

/*
Update the father and sons of a node when spliting nodes

Args:
    old_node [Node*]: [the old node]
    new_node [Node*]: [the new node]
*/
void TrapezoidalTree::UpdateFatherSon(Node* old_node, Node* new_node)
{
    if(old_node->fathers.size() == 0)
    {
        this->root = new_node;
    }
    else 
    {
        for(int i = 0; i < old_node->fathers.size(); i ++)
        {
            Node* father = old_node->fathers[i];
            if(father->left_son == old_node)
            {
                father->left_son = new_node;
            }
            if(father->right_son == old_node)
            {
                father->right_son = new_node;
            }
            new_node->fathers.push_back(father);
        }
    }
}

/*
Judge that whether two trap nodes can merge

Args:
    left [TrapNode*]: [the left trapnode]
    right [TrapNode*]: [the right trapnode]

Returns:
    result [bool]: [can merge or not]
*/
bool TrapezoidalTree::JudgeTrapNodeMerge(TrapNode* left, TrapNode* right)
{
    bool result = 0;
    if(left->trap->up->id == right->trap->up->id && left->trap->down->id == right->trap->down->id)
    {
        result = 1;
    }
    return result;
}

/*
Merge the adjacent trap nodes
*/
void TrapezoidalTree::MergeTrapNode(TrapNode* left, TrapNode* right)
{
    //update left trap
    left->trap->right = right->trap->right;

    //update father and son
    this->UpdateFatherSon(right, left);

    //update topology
    left->right_down = right->right_down;
    this->UpdateNeighborTopology(right, left, right->right_down);
    left->right_up = right->right_up;
    this->UpdateNeighborTopology(right, left, right->right_up);

    //delete right
    delete(right);
}

/*
Insert a new line into the trapezoidal tree

Args:
    trap [TrapNode*]: [the trap containing the line's left point]
    l [Line*]: [the new line to be inserted]
*/
void TrapezoidalTree::InsertLine(TrapNode* trap, Line* l)
{
    //find all relavant trapezoidal nodes
    vector<TrapNode*> split_nodes;
    split_nodes.clear();
    while(trap != NULL)
    {
        split_nodes.push_back(trap);
        //judge whether the line will cross the trap
        double line_right = l->t->x;
        double trap_right = trap->trap->right;
        if(line_right < trap_right) //not cross
        {
            break;
        }
        else //cross
        {
            //judge the next cross: right up or right down
            TrapNode* right_down_node = trap->right_down;
            if(right_down_node == NULL)
            {
                trap = trap->right_up;
            }
            else 
            {
                double y_line = l->GetIntersectionY(trap_right);
                double y_right_down = trap->right_down->trap->up->GetIntersectionY(trap_right);
                if(y_line >= y_right_down)
                {
                    trap = trap->right_up;
                }
                else 
                {
                    trap = trap->right_down;
                }
            }
        }
    }


    //split the nodes
    //if only one node: split into four pieces
    if(split_nodes.size() == 1)
    {
        //split the traps
        TrapNode* old_trap_node = split_nodes[0];
        Trap* old_trap = old_trap_node->trap;
        double left_x = old_trap->left;
        double right_x = old_trap->right;
        double split_left_x = l->s->x;
        double split_right_x = l->t->x;
        Line* up_line = old_trap->up;
        Line* down_line = old_trap->down;
        Trap* L_trap = new Trap(left_x, split_left_x, up_line, down_line);
        Trap* A_trap = new Trap(split_left_x, split_right_x, up_line, l);
        Trap* B_trap = new Trap(split_left_x, split_right_x, l, down_line);
        Trap* R_trap = new Trap(split_right_x, right_x, up_line, down_line);
            
        //build nodes
        PointNode* p = new PointNode(l->s);
        PointNode* q = new PointNode(l->t);
        LineNode* r = new LineNode(l);
        TrapNode* L = new TrapNode(L_trap);
        TrapNode* A = new TrapNode(A_trap);
        TrapNode* B = new TrapNode(B_trap);
        TrapNode* R = new TrapNode(R_trap);

        //update father and son
        this->UpdateFatherSon(old_trap_node, p);
        p->left_son = L;
        p->right_son = q;
        L->fathers.push_back(p);
        q->fathers.push_back(p);
        q->left_son = r;
        q->right_son = R;
        r->fathers.push_back(q);
        R->fathers.push_back(q);
        r->left_son = A;
        r->right_son = B; 
        A->fathers.push_back(r);
        B->fathers.push_back(r);

        //update trapnode topology
        L->left_up = old_trap_node->left_up;
        L->left_down = old_trap_node->left_down;
        L->right_up = A;
        L->right_down = B; 
        A->left_up = L;
        A->left_down = L;
        A->right_up = R;
        A->right_down = R;
        B->left_up = L;
        B->left_down = L;
        B->right_up = R;
        B->right_down = R;            
        R->left_up = A;
        R->left_down = B; 
        R->right_up = old_trap_node->right_up;
        R->right_down = old_trap_node->right_down;
        this->UpdateNeighborTopology(old_trap_node, L, old_trap_node->left_up);
        this->UpdateNeighborTopology(old_trap_node, L, old_trap_node->left_down);
        this->UpdateNeighborTopology(old_trap_node, R, old_trap_node->right_up);
        this->UpdateNeighborTopology(old_trap_node, R, old_trap_node->right_down);

        split_nodes.clear();
        delete(old_trap_node);
    }

    //many new nodes
    else 
    {
        //build the two point nodes and the left and right nodes
        //node
        TrapNode* old_left_node = split_nodes[0];
        TrapNode* old_right_node = split_nodes[split_nodes.size() - 1];
        double l_left = old_left_node->trap->left;
        double l_right = l->s->x;
        double r_left = l->t->x;
        double r_right = old_right_node->trap->right;
        Line* l_up = old_left_node->trap->up;
        Line* l_down = old_left_node->trap->down;
        Line* r_up = old_right_node->trap->up;
        Line* r_down = old_right_node->trap->down;
        Trap* L_trap = new Trap(l_left, l_right, l_up, l_down);
        Trap* R_trap = new Trap(r_left, r_right, r_up, r_down);
        PointNode* p = new PointNode(l->s);
        PointNode* q = new PointNode(l->t);
        TrapNode* L = new TrapNode(L_trap);
        TrapNode* R = new TrapNode(R_trap);
        //tree
        this->UpdateFatherSon(old_left_node, p);
        this->UpdateFatherSon(old_right_node, q);
        p->left_son = L;
        L->fathers.push_back(p);
        q->right_son = R;
        R->fathers.push_back(q);
        //topology
        L->left_down = old_left_node->left_down;
        L->left_up = old_left_node->left_up;
        R->right_down = old_right_node->right_down;
        R->right_up = old_right_node->right_up;
        this->UpdateNeighborTopology(old_left_node, L, old_left_node->left_up);
        this->UpdateNeighborTopology(old_left_node, L, old_left_node->left_down);
        this->UpdateNeighborTopology(old_right_node, R, old_right_node->right_up);
        this->UpdateNeighborTopology(old_right_node, R, old_right_node->right_down);

        //visit all the trap nodes to construct the A and B splits
        list<LineNode*> new_rs;
        list<TrapNode*> new_As;
        list<TrapNode*> new_Bs;
        new_rs.clear();
        new_As.clear();
        new_Bs.clear();
        for(int i = 0; i < split_nodes.size(); i ++)
        {
            //traps and nodes
            TrapNode* old_node = split_nodes[i];
            double left = old_node->trap->left;
            if(i == 0)
            {
                left = l_right;
            }
            double right = old_node->trap->right;
            if(i == split_nodes.size() - 1)
            {
                right = r_left;
            }
            Line* up = old_node->trap->up;
            Line* down = old_node->trap->down;
            Trap* A_trap = new Trap(left, right, up, l);
            Trap* B_trap = new Trap(left, right, l, down);
            TrapNode* A = new TrapNode(A_trap);
            TrapNode* B = new TrapNode(B_trap);
            LineNode* r = new LineNode(l);
            //tree
            if(i == 0)
            {
                p->right_son = r;
                r->fathers.push_back(p);
            }
            else if(i == split_nodes.size() - 1)
            {
                q->left_son = r;
                r->fathers.push_back(q);
            }
            else 
            {
                this->UpdateFatherSon(old_node, r);
            }
            r->left_son = A;
            r->right_son = B; 
            A->fathers.push_back(r);
            B->fathers.push_back(r);
            //store
            new_rs.push_back(r);
            new_As.push_back(A);
            new_Bs.push_back(B);
        }
        //topology, only inherit
        list<TrapNode*>::iterator ia = new_As.begin(); 
        list<TrapNode*>::iterator ib = new_Bs.begin(); 
        for(int i = 0; i < split_nodes.size(); i ++)
        {
            TrapNode* current_a = (*ia);
            TrapNode* current_b = (*ib);
            if(i != 0)
            {
                if(split_nodes[i]->left_up != split_nodes[i - 1])
                {
                    current_a->left_up = split_nodes[i]->left_up;
                    this->UpdateNeighborTopology(split_nodes[i], current_a, split_nodes[i]->left_up);
                }
                if(split_nodes[i]->left_down != split_nodes[i - 1])
                {
                    current_b->left_down = split_nodes[i]->left_down;
                    this->UpdateNeighborTopology(split_nodes[i], current_b, split_nodes[i]->left_down);
                }
            }
                
            if(i != split_nodes.size() - 1)
            {
                if(split_nodes[i]->right_up != split_nodes[i + 1])
                {
                    current_a->right_up = split_nodes[i]->right_up;
                    this->UpdateNeighborTopology(split_nodes[i], current_a, split_nodes[i]->right_up);
                }
                if(split_nodes[i]->right_down != split_nodes[i + 1])
                {
                    current_b->right_down = split_nodes[i]->right_down;
                    this->UpdateNeighborTopology(split_nodes[i], current_b, split_nodes[i]->right_down);
                }
            }
            ia ++;
            ib ++;
        }

        //merge the redundant As and Bs
        ia = new_As.begin();
        while(ia != new_As.end())
        {
            while(ia != new_As.end())
            {
                ia ++;
                list<TrapNode*>::iterator next_ia = ia;
                ia --;
                if(next_ia == new_As.end())
                {
                    break;
                }

                TrapNode* this_a = (*ia);
                TrapNode* next_a = (*next_ia);
                if(this->JudgeTrapNodeMerge(this_a, next_a))
                {
                    this->MergeTrapNode(this_a, next_a);
                    new_As.erase(next_ia);
                }
                else 
                {
                    break;
                }
            }
            ia ++;
        }

        ib = new_Bs.begin();
        while(ib != new_Bs.end())
        {
            while(ib != new_Bs.end())
            {
                ib ++;
                list<TrapNode*>::iterator next_ib = ib;
                ib --;
                if(next_ib == new_Bs.end())
                {
                    break;
                }

                TrapNode* this_b = (*ib);
                TrapNode* next_b = (*next_ib);
                if(this->JudgeTrapNodeMerge(this_b, next_b))
                {
                    this->MergeTrapNode(this_b, next_b);
                    new_Bs.erase(next_ib);
                }
                else 
                {
                    break;
                }
            }
            ib ++;
        }

        //set topology after merging
        ia = new_As.begin();
        while(ia != new_As.end())
        {
            TrapNode* current_node = (*ia);
            if(ia == new_As.begin())
            {
                L->right_up = current_node;
                current_node->left_up = L;
                current_node->left_down = L;
            }
            else
            {
                ia --;
                TrapNode* last_node = (*ia);
                ia ++;
                last_node->right_down = current_node;
                current_node->left_down = last_node;
                if(last_node->right_up == NULL)
                {
                    last_node->right_up = current_node;
                }
                if(current_node->left_up == NULL)
                {
                    current_node->left_up = last_node;
                }
            }
            ia ++;
        }
        ia = new_As.end();
        ia --;
        (*ia)->right_up = R;
        (*ia)->right_down = R;
        R->left_up = (*ia);

        ib = new_Bs.begin();
        while(ib != new_Bs.end())
        {
            TrapNode* current_node = (*ib);
            if(ib == new_Bs.begin())
            {
                L->right_down = current_node;
                current_node->left_up = L;
                current_node->left_down = L;
            }
            else
            {
                ib --;
                TrapNode* last_node = (*ib);
                ib ++;
                last_node->right_up = current_node;
                current_node->left_up = last_node;
                if(last_node->right_down == NULL)
                {
                    last_node->right_down = current_node;
                }
                if(current_node->left_down == NULL)
                {
                    current_node->left_down = last_node;
                }
            }
            ib ++;
        }
        ib = new_Bs.end();
        ib --;
        (*ib)->right_up = R;
        (*ib)->right_down = R;
        R->left_down = (*ib);


        for(int i = 0; i < split_nodes.size(); i ++)
        {
            delete(split_nodes[i]);
        }
        split_nodes.clear();
        new_As.clear();
        new_Bs.clear();
        new_rs.clear();
    }
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
    (*best_dist) = 2000000.0;
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
    while(current_place >= 0 && current_place < lines_sx.size() && lines_sx[current_place]->s->x == target_x)
    {
        
        double current_y = lines_sx[current_place]->s->y;
        
        double dist_y = current_y - target_y;

        //judge vertical lines
        if(lines_sx[current_place]->vertical == 1)
        {
            if(target_y >= lines_sx[current_place]->s->y && target_y <= lines_sx[current_place]->t->y)
            {
                dist_y = 0;
            }
        }

        if(dist_y >= 0 && dist_y < *best_dist)
        {
            *best_dist = dist_y;
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
    current_place = start;
    while(current_place >= 0 && current_place < lines_tx.size() && lines_tx[current_place]->t->x == target_x)
    {
        double current_y = lines_tx[current_place]->t->y;
        double dist_y = current_y - target_y;

        //judge vertical lines
        if(lines_tx[current_place]->vertical == 1)
        {
            if(target_y >= lines_tx[current_place]->s->y && target_y <= lines_tx[current_place]->t->y)
            {
                dist_y = 0;
            }
        }

        if(dist_y >= 0 && dist_y < *best_dist)
        {
            *best_dist = dist_y;
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
    all_lines_sx.clear();
    all_lines_tx.clear();
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

    //randomly shuffle the lines
    mt19937 random_generator;
    random_generator.seed(14530529); //set seed to reproduce
    shuffle(normal_lines.begin(), normal_lines.end(), random_generator);

    //main
    TrapezoidalTree* the_tree = new TrapezoidalTree(normal_lines);
    the_tree->BuildTree();
    for(int i = 0; i < query_points.size(); i ++)
    {
        TrapNode* best_trap_tree = the_tree->SearchPoint(the_tree->root, query_points[i]);
        Line* best_line_tree = best_trap_tree->trap->up;
        double best_dist_tree = best_line_tree->GetIntersectionY(query_points[i]->x) - query_points[i]->y;

        double best_dist_special_judge;
        Line* best_line_special_judge = SpecialJudge(query_points[i], all_lines_sx, all_lines_tx, &best_dist_special_judge);
        Line* best_line;
        if(best_dist_tree < best_dist_special_judge)
        {
            best_line = best_line_tree;
        }
        else 
        {
            best_line = best_line_special_judge;
        }
        long long best_id = best_line->id;
        if(best_id >= 1 && best_id <= n)
        {
            printf("%lld\n", best_id);
        }
        else 
        {
            printf("N\n");
        }
    }
    return 0;
}