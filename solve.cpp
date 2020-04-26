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

using my_pair_t = std::pair<point, size_t>;
using my_container_t = std::vector<my_pair_t>;

class MyHashFunction
{
public:
    size_t operator()(const point &p) const
    {
        return (hash<int>()(p.first)) ^ (hash<int>()(p.second) << 1);
    }
};

class MyEqualFunction
{
public:
    bool operator()(const point &p1, const point &p2) const
    {
        return ((p1.first == p2.first) && (p1.second == p2.second));
    }
};

bool getUnvisitedNeighbors(point &neighbor,
                           unordered_set<point, MyHashFunction, MyEqualFunction> &visited,
                           point currPosition,
                           Maze &m)
{
    for (uint16_t i = 0; i < 4; i++)
    {
        neighbor = currPosition + moveIn(i);
        if (m.can_go(i, currPosition.first, currPosition.second) &&
            visited.find(neighbor) == visited.end())
        {
            return true;
        }
    }
    return false;
}

path getShortestpath(Maze &m, point start, point end)
{
    list<point> path;
    point newPosition;
    size_t cost;
    my_pair_t currPosition;
    unordered_map<point, int, MyHashFunction, MyEqualFunction> cellCost;

    for (int i = 0; i < m.rows() ; i++)
    {
        for (int j = 0; j < m.columns(); j++)
        {
            point p = make_pair(i, j);
            cellCost[p] = INT_MAX;
        }
    }

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;
    unordered_map<point, point, MyHashFunction, MyEqualFunction> prev;

    auto my_comp = [](const my_pair_t &e1, const my_pair_t &e2) { return e1.second > e2.second; };

    std::priority_queue<my_pair_t,
                        my_container_t,
                        decltype(my_comp)> queue(my_comp);

    queue.emplace(start, 0);

    while (!queue.empty())
    {
        currPosition = queue.top();
        point curr = queue.top().first;
        queue.pop();

        if (curr == end)
        {
            break;
        }

        if (visited.find(curr) == visited.end())
        {
            visited.emplace(curr);
        }

        for (int16_t i = 0; i < 4; i++)
        {
            if (m.can_go(i, curr.first, curr.second) && visited.find(curr + moveIn(i)) == visited.end())
            {
                newPosition = curr + moveIn(i);
                cost = currPosition.second + m.cost(curr, i);
                if (cellCost[newPosition] > cost)
                {
                    cellCost[newPosition] = cost;
                    queue.emplace(newPosition, cost);
                    prev[newPosition] = curr;
                }
            }
        }
    }

    point currPt = end;

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

path solve_dfs(Maze &m, int rows, int cols)
{
    stack<point> stacklist;
    list<point> path;
    point currPosition;
    
    point start = make_pair(0, 0);
    point end = make_pair(rows - 1, cols - 1);

    stacklist.push(start);
    path.push_back(start);

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;
    point newPosition;
    bool found = false;

    while (true)
    {
        currPosition = stacklist.top();
        if (currPosition == end)
            break;

        if (visited.find(currPosition) == visited.end())
            visited.emplace(currPosition);

        found = getUnvisitedNeighbors(newPosition, visited, currPosition, m);

        if (found)
        {
            stacklist.push(newPosition);
            path.push_back(newPosition);
        }
        else
        {
            stacklist.pop();
            path.pop_back();
        }
    }
    return path;
}

path solve_bfs(Maze &m, int rows, int cols)
{
    queue<point> queuelist;
    point start = make_pair(0, 0);
    point end = make_pair(rows - 1, cols - 1);

    list<point> path;
    point currPosition;
    point newPosition;

    queuelist.emplace(start);

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;
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
            newPosition = currPosition + moveIn(i);
            if (m.can_go(i, currPosition.first, currPosition.second) && visited.find(newPosition) == visited.end())
            {
                visited.emplace(newPosition);
                queuelist.emplace(newPosition);
                prev[newPosition] = currPosition;
            }
        }
    }

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

path solve_dijkstra(Maze &m, int rows, int cols)
{
    point start = make_pair(0, 0);
    point end = make_pair(rows - 1, cols - 1);

    list<point> path = getShortestpath(m, start, end);

    return path;
}

path solve_tour(Maze &m, int rows, int cols)
{
    point start = make_pair(rows/2, cols/2);

    point corner1 = make_pair(0, 0);
    point corner2 = make_pair(0, cols - 1);
    point corner3 = make_pair(rows - 1, cols - 1);
    point corner4 = make_pair(rows - 1, 0);

    list<point> path = getShortestpath(m, start, corner1);

    list<point> path1 = getShortestpath(m, corner1, corner2);
    path1.pop_front();

    list<point> path2 = getShortestpath(m, corner2, corner3);
    path2.pop_front();

    list<point> path3 = getShortestpath(m, corner3, corner4);
    path3.pop_front();

    list<point> path4 = getShortestpath(m, corner4, start);
    path4.pop_front();

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
