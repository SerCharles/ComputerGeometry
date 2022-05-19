#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

#define SPLIT_THRESHOLD 50 //< this threshold: leaf node

/*
Definition of a point
Args:
    x [double]: the x coordinate of the point
    y [double]: the y coordinate of the point
*/
struct Point
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


//used in getting the median
bool CompareX(Point a, Point b)
{
    return a.x < b.x;
}
bool CompareY(Point a, Point b)
{
    return a.y < b.y;
}

class BoundingBox
{
public:
    int left;
    int right;
    int down;
    int up;

    BoundingBox()
    {
        left = 0;
        right = 0;
        down = 0;
        up = 0;
    }

    BoundingBox(int left_, int right_, int down_, int up_)
    {
        left = left_;
        right = right_;
        down = down_;
        up = up_;
    }

    BoundingBox(vector<Point>& points)
    {
        left = 2147483647;
        right = -2147483648;
        down = 2147483647;
        up = -2147483648;
        for(int i = 0; i < points.size(); i ++)
        {
            int x = points[i].x;
            int y = points[i].y;
            if(x < left)
            {
                left = x;
            }
            if(x > right)
            {
                right = x;
            }
            if(y < down)
            {
                down = y;
            }
            if(y > up)
            {
                up = y;
            }
        }
    }

    /*
    Judge whether a point p is inside the bounding box

    Args:
        p [Point]: [the point to be judged]

    Returns:
        result [bool]: [inside: 1, else: 0]
    */
    bool JudgeInside(Point p)
    {
        bool result = 0;
        if(p.x >= left && p.x <= right && p.y >= down && p.y <= up)
        {
            result = 1;
        }
        return result;
    }

    /*
    Judge whether a bounding box b is inside the bounding box

    Args:
        b [BoundingBox]: [the bounding box to be judged]

    Returns:
        result [bool]: [a contains b: 1, else: 0]
    */
    bool JudgeContain(BoundingBox b)
    {
        bool result = 0;
        if(b.left >= left && b.right <= right && b.down >= down && b.up <= up)
        {
            result = 1;
        }
        return result;
    }

    /*
    Judge whether a bounding box b intersects with the bounding box

    Args:
        b [BoundingBox]: [the bounding box to be judged]

    Returns:
        result [bool]: [a intersects with b: 1, else: 0]
    */
    bool JudgeIntersection(BoundingBox b)
    {
        bool result = 0;
        int max_left = max(left, b.left);
        int min_right = min(right, b.right);
        int max_down = max(down, b.down);
        int min_up = min(up, b.up);
        if(max_left <= min_right && max_down <= min_up)
        {
            result = 1;
        }
        return result;
    }
};

class RangeNodeX;
class RangeNodeY
{
public:
    int depth = 0; //the depth of current node
    int start = 0; //the starting id of the points in the range
    int point_number = 0; //the number of points in the range
    int y = 0; //the y value of the node 
    bool leaf = 0; //whether it is leaf node or not
    RangeNodeX* father = NULL; //the father X tree
    RangeNodeY* left_son = NULL;
    RangeNodeY* right_son = NULL;

    ~RangeNodeY()
    {
        delete(this->left_son);
        delete(this->right_son);
    }

    RangeNodeY(int depth_, RangeNodeX* father_, int start_, int point_number_);
};

class RangeNodeX
{
public:
    int depth = 0; //the depth of current node
    int point_number = 0; //the number of points in the range
    int x = 0; //the x value of the node 
    bool leaf = 0; //whether it is leaf node or not
    vector<Point> points; //the points in the range 
    RangeNodeY* y_tree = NULL; //the y tree
    RangeNodeX* left_son = NULL;
    RangeNodeX* right_son = NULL;

    ~RangeNodeX()
    {
        this->points.clear();
        delete(this->y_tree);
        delete(this->left_son);
        delete(this->right_son);
    }
    RangeNodeX(int depth_, vector<Point>& points_);
};

RangeNodeY::RangeNodeY(int depth_, RangeNodeX* father_, int start_, int point_number_)
{
    //set values
    this->depth = depth_;
    this->father = father_;
    this->start = start_;
    this->point_number = point_number_;

    //get median and split the list
    int median_index = this->start + this->point_number / 2;
    nth_element(this->father->points.begin() + this->start, this->father->points.begin() + median_index, 
        this->father->points.begin() + this->start + this->point_number, CompareY);
    this->y = this->father->points[median_index].y;

    //leaf node 
    if(this->point_number <= SPLIT_THRESHOLD)
    {
        this->leaf = 1;
        this->left_son = NULL;
        this->right_son = NULL;
    }
    //build son
    else 
    {
        this->leaf = 0;
        int left_start = this->start;
        int left_length = median_index - this->start;
        int right_start = this->start + left_length;
        int right_length = this->point_number - left_length;
        this->left_son = new RangeNodeY(this->depth + 1, this->father, left_start, left_length);
        this->right_son = new RangeNodeY(this->depth + 1, this->father, right_start, right_length);
    }
}

RangeNodeX::RangeNodeX(int depth_, vector<Point>& points_)
{
    //set values
    this->depth = depth_;
    this->point_number = points_.size();
    this->points.clear();
    for(int i = 0; i < this->point_number; i ++)
    {
        this->points.push_back(points_[i]);
    }

    //get median and split the list
    int median_index = this->point_number / 2;
    nth_element(this->points.begin(), this->points.begin() + median_index, this->points.end(), CompareX);
    this->x = this->points[median_index].x;

    //leaf node 
    if(this->point_number <= SPLIT_THRESHOLD)
    {
        this->leaf = 1;
        this->left_son = NULL;
        this->right_son = NULL;
        this->y_tree = NULL;
    }
    //build son
    else 
    {
        this->leaf = 0;
        vector<Point> left_list;
        vector<Point> right_list;
        left_list.clear();
        right_list.clear();
        for(int i = 0; i < median_index; i ++)
        {
            left_list.push_back(this->points[i]);
        }

        for(int i = median_index; i < this->point_number; i ++)
        {
            right_list.push_back(this->points[i]);
        }
        this->left_son = new RangeNodeX(this->depth + 1, left_list);
        this->right_son = new RangeNodeX(this->depth + 1, right_list);
        left_list.clear();
        right_list.clear();

        //build Y
        this->y_tree = new RangeNodeY(1, this, 0, this->point_number);
    }  
}


/*
Count the points inside the bounding box

Args:
    box [BoundingBox]: [the bounding box to be searched]
    root [RangeNodeX]: [the root point of the Y range tree]

Return:
    result [int]: [the number of total points]
*/
int CountPoints(BoundingBox& box, RangeNodeY* root)
{
    int min_y = box.down;
    int max_y = box.up;
    vector<RangeNodeY*> route_left;
    vector<RangeNodeY*> route_right;
    vector<RangeNodeY*> canonical_nodes;
    route_left.clear();
    route_right.clear();
    canonical_nodes.clear();

    //get the left route
    RangeNodeY* current_left = root;
    route_left.push_back(current_left);
    while(current_left->leaf == 0)
    {
        if(min_y <= current_left->y)
        {
            current_left = current_left->left_son;
            route_left.push_back(current_left);
        }
        else 
        {
            current_left = current_left->right_son;
            route_left.push_back(current_left);
        }
    }

    //get the right route
    RangeNodeY* current_right = root;
    route_right.push_back(current_right);
    while(current_right->leaf == 0)
    {
        if(max_y >= current_right->y)
        {
            current_right = current_right->right_son;
            route_right.push_back(current_right);
        }
        else 
        {
            current_right = current_right->left_son;
            route_right.push_back(current_right);
        }
    }

    //get lca
    int lca = 0;
    for(int i = 0; i < min(route_left.size(), route_right.size()); i ++)
    {
        if(route_left[i] == route_right[i])
        {
            lca = i;
        }
        else 
        {
            break;
        }
    }

    //get canonical nodes
    if(route_left[lca]->leaf)
    {
        canonical_nodes.push_back(route_left[lca]);
    }
    else 
    {
        for(int i = lca + 1; i < route_left.size(); i ++)
        {
            if(route_left[i]->leaf)
            {
                canonical_nodes.push_back(route_left[i]);
            }
            else if(route_left[i]->y >= min_y)
            {
                canonical_nodes.push_back(route_left[i]->right_son);
            }
        }
        for(int i = lca + 1; i < route_right.size(); i ++)
        {
            if(route_right[i]->leaf)
            {
                canonical_nodes.push_back(route_right[i]);
            }
            else if(route_right[i]->y <= max_y)
            {
                canonical_nodes.push_back(route_right[i]->left_son);
            }
        }
    }

    //count
    int result = 0;
    for(int i = 0; i < canonical_nodes.size(); i ++)
    {
        int the_result = 0;
        if(canonical_nodes[i]->leaf)
        {
            for(int j = canonical_nodes[i]->start; j < canonical_nodes[i]->start + canonical_nodes[i]->point_number; j ++)
            {
                if(box.JudgeInside(canonical_nodes[i]->father->points[j]))
                {
                    the_result += 1;
                }
            }
        }
        else 
        {
            the_result = canonical_nodes[i]->point_number;
        }
        result = result + the_result;
    }
    return result;
}

/*
Count the points inside the bounding box

Args:
    box [BoundingBox]: [the bounding box to be searched]
    root [RangeNodeX]: [the root point of the X range tree]

Return:
    result [int]: [the number of total points]
*/
int CountPoints(BoundingBox& box, RangeNodeX* root)
{
    int min_x = box.left;
    int max_x = box.right;
    vector<RangeNodeX*> route_left;
    vector<RangeNodeX*> route_right;
    vector<RangeNodeX*> canonical_nodes;
    route_left.clear();
    route_right.clear();
    canonical_nodes.clear();

    //get the left route
    RangeNodeX* current_left = root;
    route_left.push_back(current_left);
    while(current_left->leaf == 0)
    {
        if(min_x <= current_left->x)
        {
            current_left = current_left->left_son;
            route_left.push_back(current_left);
        }
        else 
        {
            current_left = current_left->right_son;
            route_left.push_back(current_left);
        }
    }

    //get the right route
    RangeNodeX* current_right = root;
    route_right.push_back(current_right);
    while(current_right->leaf == 0)
    {
        if(max_x >= current_right->x)
        {
            current_right = current_right->right_son;
            route_right.push_back(current_right);
        }
        else 
        {
            current_right = current_right->left_son;
            route_right.push_back(current_right);
        }
    }

    //get lca
    int lca = 0;
    for(int i = 0; i < min(route_left.size(), route_right.size()); i ++)
    {
        if(route_left[i] == route_right[i])
        {
            lca = i;
        }
        else 
        {
            break;
        }
    }

    //get canonical nodes
    if(route_left[lca]->leaf)
    {
        canonical_nodes.push_back(route_left[lca]);
    }
    else 
    {
        for(int i = lca + 1; i < route_left.size(); i ++)
        {
            if(route_left[i]->leaf)
            {
                canonical_nodes.push_back(route_left[i]);
            }
            else if(route_left[i]->x >= min_x)
            {
                canonical_nodes.push_back(route_left[i]->right_son);
            }
        }
        for(int i = lca + 1; i < route_right.size(); i ++)
        {
            if(route_right[i]->leaf)
            {
                canonical_nodes.push_back(route_right[i]);
            }
            else if(route_right[i]->x <= max_x)
            {
                canonical_nodes.push_back(route_right[i]->left_son);
            }
        }
    }
    
    //count
    int result = 0;
    for(int i = 0; i < canonical_nodes.size(); i ++)
    {
        int the_result = 0;
        if(canonical_nodes[i]->leaf)
        {
            for(int j = 0; j < canonical_nodes[i]->point_number; j ++)
            {
                if(box.JudgeInside(canonical_nodes[i]->points[j]))
                {
                    the_result += 1;
                }
            }
        }
        else 
        {
            the_result = CountPoints(box, canonical_nodes[i]->y_tree);
        }
        result = result + the_result;
    }
    return result;
}

int main()
{
    vector<Point> points;
    vector<BoundingBox> rectangles;
    points.clear();
    rectangles.clear();
    int n, m;
    scanf("%d %d", &n, &m);
    for(int i = 0; i < n; i ++)
    {
        int x, y;
        scanf("%d %d", &x, &y);
        Point new_point = Point(x, y);
        points.push_back(new_point);
    }
    for(int i = 0; i < m; i ++)
    {
        int left, down, right, up;
        scanf("%d %d %d %d", &left, &down, &right, &up);
        BoundingBox new_rectangle = BoundingBox(left, right, down, up);
        rectangles.push_back(new_rectangle);
    }
    RangeNodeX* range_tree = new RangeNodeX(1, points);

    for(int i = 0; i < m; i ++)
    {
        int point_number = CountPoints(rectangles[i], range_tree);
        printf("%d\n", point_number);
    }
    delete(range_tree);
    return 0;
}