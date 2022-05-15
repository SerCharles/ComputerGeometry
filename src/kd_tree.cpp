#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
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
vector<Point> points;

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

class KDTree
{
public: 
    int depth; //the depth of the current node. depth_ % 2 == 1: split x, or split y
    int start_id; //the start id of the points contained in the KDNode of all points
    int point_number; //the point numbers
    bool leaf; //whether the node is leaf or not
    BoundingBox bounding_box; //the bounding box of the node
    bool split_x;
    int threshold = 50; //if number of points <= threshold, then set it to the leaf node
    KDTree* left_son = NULL;
    KDTree* right_son = NULL;
    KDTree()
    {
        depth = 0;
        start_id = 0;
        point_number = 0;
        leaf = 0;
        left_son = NULL;
        right_son = NULL;
        split_x = 0;
    }

    //used in getting the median
    struct Compare
    {
        bool split_x; //split x or y
        bool operator() (Point a, Point b) const
        {
            if (split_x)
            {
                return a.x < b.x;
            }
            else
            {
                return a.y < b.y;
            }
        }
        Compare(bool x)
        {
            split_x = x;
        }
    };

    KDTree(int depth_, int start_id_, int point_number_, BoundingBox bounding_box_)
    {
        //init basic information
        depth = depth_;
        start_id = start_id_;
        point_number = point_number_;
        bounding_box = bounding_box_;
        split_x = (depth_ % 2 == 1) ? 1 : 0;

        if(point_number <= threshold)
        {
            leaf = 1;
            left_son = NULL;
            right_son = NULL;
        }
        else 
        {
            leaf = 0;
            int median_index = start_id + point_number / 2;
            nth_element(points.begin() + start_id, points.begin() + median_index, points.begin() + start_id + point_number, Compare(split_x));
            
            //get median
            int median = split_x ? points[median_index].x : points[median_index].y;

            //get new bounding box
            BoundingBox bounding_box_left = bounding_box;
            BoundingBox bounding_box_right = bounding_box;
            if(split_x)
            {
                bounding_box_left.right = median;
                bounding_box_right.left = median;
            }
            else 
            {
                bounding_box_left.up = median;
                bounding_box_right.down = median;
            }

            //get sons 
            left_son = new KDTree(depth + 1, start_id, median_index - start_id, bounding_box_left);
            right_son = new KDTree(depth + 1, median_index, start_id + point_number - median_index, bounding_box_right);
        }
    }

    /*
    Get the number of points inside b

    Args:
        b [BoundingBox]: [the bounding box to be queried]

    Returns:
        result [int]: [the number of points inside b]
    */
    int GetPointNumber(BoundingBox& b)
    {
        int result = 0;
        //current tree is inside the query box
        if(b.JudgeContain(bounding_box))
        {
            result = point_number;
            return result;
        }
        //not intersect
        else if(b.JudgeIntersection(bounding_box) == 0)
        {
            result = 0;
            return result;
        }
        //intersect but not contain
        else 
        {
            //leaf: brute force
            if(leaf)
            {
                for(int i = start_id; i < start_id + point_number; i ++)
                {
                    if(b.JudgeInside(points[i]))
                    {
                        result += 1;
                    }
                }
                return result;
            }
            else 
            {
                result = left_son->GetPointNumber(b) + right_son->GetPointNumber(b);
                return result;
            }
        }
        return result;
    }   

};



int main()
{
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

    BoundingBox global_bounding_box = BoundingBox(points);
    KDTree* my_tree = new KDTree(1, 0, n, global_bounding_box);

    for(int i = 0; i < m; i ++)
    {
        int point_number = my_tree->GetPointNumber(rectangles[i]);
        printf("%d\n", point_number);
    }
    return 0;
}