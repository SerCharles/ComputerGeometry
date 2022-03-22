#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "point.hpp"
#include "gs.hpp"
using namespace std;



int main()
{
    //input
    int total_number = 0;
    int big_number = 1000000007;
    vector<Point> total_points;
    total_points.clear();
    scanf("%d", &total_number);
    for(int i = 0; i < total_number; i ++)
    {
        int x = 0, y = 0;
        scanf("%d %d", &x, &y);
        Point new_point(x, y);
        total_points.push_back(new_point);
    }
    int result = GrahamScan(total_points);
    printf("%d", result);
    return 0;
}