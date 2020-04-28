#include "maze.h"
#include "path.h"
#include <queue>
#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include <iostream>
#include <climits>
#include <sstream>
#include <stack>
#include <unordered_set>
#include <unordered_map>
using namespace std;

path solve_dfs(Maze &m, int rows, int cols);
path solve_bfs(Maze &m, int rows, int cols);
path solve_dijkstra(Maze &m, int rows, int cols);
path solve_tour(Maze &m, int rows, int cols);

/// pair of point and integer
using my_pair_t = std::pair<point, size_t>;

/// containter to hold vector of pairs
using my_container_t = std::vector<my_pair_t>;


/* Structure definition for hash function
    to create a hash value for a point(pair of int)
 */
struct MyHashFunction
{
    size_t operator()(const point &p) const
    {
        return (hash<int>()(p.first)) ^ (hash<int>()(p.second) << 1);
    }
};

/* Structure definition for equal function
    to compare 2 points(pair of int) 
    return true if first element of point1 is equal to first element of point2 
    and second element of point1 is equal to secons element of point2 
*/
struct MyEqualFunction
{
    bool operator()(const point &p1, const point &p2) const
    {
        return ((p1.first == p2.first) && (p1.second == p2.second));
    }
};


/** 
 * check for unvisted Neighbors from the current location and return true if there is atleat 1 such neighbor exits.
 * @param m                 input maze
 * @param newNeighbor       neighboring cell of the current cell
 * @param visited           unordered set of points to hold the points which are alreday visited
 * @param currPosition      current location from where we are trying to find it's Neighbors 
 * 
 */
bool findUnvisitedNeighbors(Maze &m,
                           point &newNeighbor,
                           unordered_set<point, MyHashFunction, MyEqualFunction> &visited,
                           point currPosition)
{
    // checking for i = 0 to i < 4 since we can move only in 4 directions(up, right, left and down)
    for (uint16_t i = 0; i < 4; i++)
    {
        // gets the next possible location based on the value of i (up, right, left and down)
        newNeighbor = currPosition + moveIn(i);
        if (m.can_go(i, currPosition.first, currPosition.second) &&
            visited.find(newNeighbor) == visited.end())
        {
            // true if there exits a neighbor and not visited
            return true;
        }
    }
    return false;
}

/**
 * gets the shortest path from the given start point to given end point in a given maze.
 * retuns the shortest path i.e list of points which joins the start point to end points
 * @param m         input maze
 * @param start     strating point from where path finding starts
 * @param end       end point where path finding loop stops
 */
path getShortestpath(Maze &m, point start, point end)
{
    list<point> path;
    point newNeighbor;

    // to store the cost of moving(difference of height of the rooms) from one point to another
    size_t cost;
    my_pair_t currPosition;

    // A map to have mapping between cell and cost(height of cell)
    unordered_map<point, int, MyHashFunction, MyEqualFunction> cellCost;

    // initially storing the max value as cost(height) for each cell,
    // later compare with the max value and the actual cost(height) and store the lesser vaue
    for (int i = 0; i < m.rows() ; i++)
    {
        for (int j = 0; j < m.columns(); j++)
        {
            point p = make_pair(i, j);
            cellCost[p] = INT_MAX;
        }
    }

    //hash set to store the non visited non neighbors of a cell
    unordered_set<point, MyHashFunction, MyEqualFunction> visited;

    //hash map to have mapping between the current cell and it's previous cell, 
    // which will help us in finding the path from start to end point
    unordered_map<point, point, MyHashFunction, MyEqualFunction> prev;

    // lambda function to to compare cost(height) of 2 point and helps in storing the pair in assending order in priority queue
    auto my_comp = [](const my_pair_t &e1, const my_pair_t &e2) { return e1.second > e2.second; };

    // A priority Queue with custome type, holds a pair of point and cost value,
    // with custome comparator which compares the cost of 2 points 
    std::priority_queue<my_pair_t,
                        my_container_t,
                        decltype(my_comp)> queue(my_comp);

    // loading the queue with starting point which has cost(height) is zero
    queue.emplace(start, 0);

    while (!queue.empty())
    {
        currPosition = queue.top();
        point curr = queue.top().first;
        queue.pop();

        // when current location is queue to end means we have reached the end point so stop the execution.
        if (curr == end)
        {
            break;
        }

        // if we don't find the current point in visited set, mark it as visited by adding it to visited set.
        if (visited.find(curr) == visited.end())
        {
            visited.emplace(curr);
        }

        // loop though 4 times since there are 4 possible movements(up, right, left, down)
        for (int16_t i = 0; i < 4; i++)
        {
            // if there is a possible move and if it's not visited expore the Neighbor
            // we need to check for unvisted cell to avoid infinate loops
            if (m.can_go(i, curr.first, curr.second) && visited.find(curr + moveIn(i)) == visited.end())
            {   
                // new cell point based on the direction of movement
                newNeighbor = curr + moveIn(i);

                // adding the current cell cost with new cell(Neighbor) cost
                // so we know the cost of moving from current cell to new cell
                cost = currPosition.second + m.cost(curr, i);

                // if the cost of new cell is more then we replace with new cost,
                // which means the new cost has the shorted distance
                if (cellCost[newNeighbor] > cost)
                {
                    cellCost[newNeighbor] = cost;
                    queue.emplace(newNeighbor, cost);
                    prev[newNeighbor] = curr;
                }
            }
        }
    }

    point currPt = end;

    // after we expore the possible movements and cost of moving from start to end points,
    // we need to join the points which connects start and end point through backtracking 
    // with the help of mapping information between current and previsous point stored in prev map
    // and store the information in path , since we are doing backtracking and storing points from end point,
    // we are pushing points from front to get the points in order
    while (prev.find(currPt) != prev.end() || currPt != start)
    {
        point prevPT = prev[currPt];
        path.push_front(currPt);
        currPt = prevPT;
    }
    path.push_front(start);
    return path;
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cerr << "usage:\n"
             << "./maze option rows cols\n"
             << " options:\n"
             << "  -dfs: depth first search (backtracking)\n"
             << "  -bfs: breadth first search\n"
             << "  -dij: dijkstra's algorithm\n"
             << "  -tour: all corners tour\n"
             << "  -basic: run dfs, bfs, and dij\n"
             << "  -advanced: run dfs, bfs, dij and tour" << endl;
        return 0;
    }
    string opt(argv[1]);

    int rows, cols;
    stringstream s;
    s << argv[2] << " " << argv[3];
    s >> rows >> cols;

    // construct a new random maze;
    Maze m(rows, cols);

    // print the initial maze out
    cout << "Initial maze" << endl;
    m.print_maze(cout, opt == "-dij" || opt == "-tour");

    if (opt == "-dfs")
    {
        cout << "\nSolved dfs" << endl;
        path p = solve_dfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);
    }

    if (opt == "-bfs")
    {
        cout << "\nSolved bfs" << endl;
        path p = solve_bfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);
    }

    if (opt == "-dij")
    {
        cout << "\nSolved dijkstra" << endl;
        path p = solve_dijkstra(m, rows, cols);
        m.print_maze_with_path(cout, p, true, false);
    }

    if (opt == "-tour")
    {
        cout << "\nSolved all courners tour" << endl;
        path p = solve_tour(m, rows, cols);
        m.print_maze_with_path(cout, p, true, true);
    }
    if (opt == "-basic")
    {
        cout << "\nSolved dfs" << endl;
        path p = solve_dfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);

        cout << "\nSolved bfs" << endl;
        p = solve_bfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);

        cout << "\nSolved dijkstra" << endl;
        p = solve_dijkstra(m, rows, cols);
        m.print_maze_with_path(cout, p, true, false);
    }
    if (opt == "-advanced")
    {
        cout << "\nSolved dfs" << endl;
        path p = solve_dfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);

        cout << "\nSolved bfs" << endl;
        p = solve_bfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);

        cout << "\nSolved dijkstra" << endl;
        p = solve_dijkstra(m, rows, cols);
        m.print_maze_with_path(cout, p, true, false);

        cout << "\nSolved all courners tour" << endl;
        p = solve_tour(m, rows, cols);
        m.print_maze_with_path(cout, p, true, true);
    }
}

/**
 * find the shortest path between start point(top left corner) to end point(bottom right corner) using depth first search.
 * and return the shorted path
 * @param m     input maze
 * @param rows  number of rows in a maze
 * @param cols  number of columns in a maze
 */ 
path solve_dfs(Maze &m, int rows, int cols)
{
    //to store the points and pop the latest visited element for next exploration
    stack<point> stacklist;
    list<point> path;
    point currPosition;
    
    // since always starting from point(0,0) top left corner
    point start = make_pair(0, 0);

    // always path ends at bottom right corner
    point end = make_pair(rows - 1, cols - 1);

    stacklist.push(start);
    path.push_back(start);

    //hash set to store the visited cell information
    unordered_set<point, MyHashFunction, MyEqualFunction> visited;

    point newNeighbor;
    bool found = false;

    while (true)
    {
        currPosition = stacklist.top();
        if (currPosition == end)
            break;

        // if the current cell is not visited mark it as visited,
        // we need this information to avoid visiting the visted cells and avoid infinite loops
        if (visited.find(currPosition) == visited.end())
            visited.emplace(currPosition);

        found = findUnvisitedNeighbors(m, newNeighbor, visited, currPosition);

        // if a unvisted neighbor exits push it to stack and also add to path list
        if (found)
        {
            stacklist.push(newNeighbor);
            path.push_back(newNeighbor);
        }
        else
        {
            stacklist.pop();
            path.pop_back();
        }
    }
    return path;
}

/**
 * find the shortest path between start point(top left corner) to end point(bottom right corner) using breadth first search.
 * and return the shorted path
 * @param m     input maze
 * @param rows  number of rows in a maze
 * @param cols  number of columns in a maze
 */ 
path solve_bfs(Maze &m, int rows, int cols)
{
    // to store the points and pop the first visited element for next exploration
    queue<point> queuelist;

    point start = make_pair(0, 0);
    point end = make_pair(rows - 1, cols - 1);

    list<point> path;
    point currPosition;
    point newNeighbor;

    // add starting point to queue 
    queuelist.emplace(start);

    // hash set to store the visited points so that we can avoid inifinate loop 
    unordered_set<point, MyHashFunction, MyEqualFunction> visited;

    // hash map for mapping between current point and it's previous point
    unordered_map<point, point, MyHashFunction, MyEqualFunction> prev;

    while (!queuelist.empty())
    {
        currPosition = queuelist.front();
        if (currPosition == end)
        {
            break;
        }
        queuelist.pop();

        for (uint16_t i = 0; i < 4; i++)
        {
            newNeighbor = currPosition + moveIn(i);
            if (m.can_go(i, currPosition.first, currPosition.second) && visited.find(newNeighbor) == visited.end())
            {
                visited.emplace(newNeighbor);
                queuelist.emplace(newNeighbor);
                prev[newNeighbor] = currPosition;
            }
        }
    }

    // after reaching the end point we have point and it's previous point mapping,
    // we do backtracking and connects the points which gives us the shortes path from start to end point.
    point currPt = end;
    while (prev.find(currPt) != prev.end() && currPt != start)
    {
        point prevPT = prev[currPt];
        path.push_front(currPt);
        currPt = prevPT;
    }
    path.push_front(start);
    return path;
}

/**
 * find the shortest path between start point(top left corner) to end point(bottom right corner) using dijkistra's shortest path.
 * @param m     input maze
 * @param rows  number of rows in a maze
 * @param cols  number of columns in a maze
 */
path solve_dijkstra(Maze &m, int rows, int cols)
{
    point start = make_pair(0, 0);
    point end = make_pair(rows - 1, cols - 1);

    list<point> path = getShortestpath(m, start, end);

    return path;
}

/**
 * find the shortest path between from middle of the maze to all 4 corner of the maze.
 * @param m     input maze
 * @param rows  number of rows in a maze
 * @param cols  number of columns in a maze
 */
path solve_tour(Maze &m, int rows, int cols)
{
    // starting the tour from middle of the maze
    point start = make_pair(rows/2, cols/2);

    // storing each corner points 
    point corner1 = make_pair(0, 0);
    point corner2 = make_pair(0, cols - 1);
    point corner3 = make_pair(rows - 1, cols - 1);
    point corner4 = make_pair(rows - 1, 0);

    // gets the shortes path from middle of the maze to corner1 (0,0)
    list<point> path = getShortestpath(m, start, corner1);

    // gets the shortes path from corner1 (0,0) to corner2 (0, col-1)
    list<point> path1 = getShortestpath(m, corner1, corner2);
    //since path's last point and path1's first point is same, will have duplicate point so removing 1 point
    path1.pop_front(); 

    // gets the shortes path from corner2 corner2 (0, col-1) to corner3 (row-1, col-1)
    list<point> path2 = getShortestpath(m, corner2, corner3);
    path2.pop_front(); 

    //gets the shortes path from corner2 corner3 (row-1, col-1) to corner2 (row-1, 0)
    list<point> path3 = getShortestpath(m, corner3, corner4);
    path3.pop_front(); 

    //gets the shortes path from corner2 corner3 (row-1, 0) to start (middle of the maze)
    list<point> path4 = getShortestpath(m, corner4, start);
    path4.pop_front(); 

    // combine all the paths
    while (!path1.empty())
    {
        path.push_back(path1.front());
        path1.pop_front();
    }
    while (!path2.empty())
    {
        path.push_back(path2.front());
        path2.pop_front();
    }
    while (!path3.empty())
    {
        path.push_back(path3.front());
        path3.pop_front();
    }
     while (!path4.empty())
    {
        path.push_back(path4.front());
        path4.pop_front();
    }
 
    return path;
}
