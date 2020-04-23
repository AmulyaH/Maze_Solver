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
    // id is returned as hash function
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
    point current = make_pair(0, 0);
    stacklist.push(current);
    path.push_back(current);

    int row = m.rows();
    int col = m.columns();
    point currPosition;
    point end = make_pair(row - 1, col - 1);
    bool flag = false;

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;

    while (true)
    {
        currPosition = stacklist.top();
        if (currPosition == end)
            break;
        point newPosition;

        if (visited.find(currPosition) == visited.end())
            visited.emplace(currPosition);

        bool found = getUnvisitedNeighbors(newPosition, visited, currPosition, m);

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

    return list<point>();
}

path solve_dijkstra(Maze &m, int rows, int cols)
{
    queue<point> queuelist;
    list<point> path;
    point current = make_pair(0, 0);
    //queuelist.push(current);
    path.push_back(current);
    int row = m.rows();
    int col = m.columns();
    point end = make_pair(row - 1, col - 1);
    my_pair_t currPosition;

    unordered_map<point, int, MyHashFunction, MyEqualFunction> cellCost;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            point p = make_pair(i, j);
            cellCost[p] = INT_MAX;
        }
    }

    /*
    vector<vector<int>> vec;

    for(int i=0; i< col ; i++)
    {
        vec.push_back(vector<int>());
        for(int j=0; j< row ; j++)
        {
            vec[i].push_back(INT_MAX);
        }
    }
    */

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;
    unordered_map<point,point, MyHashFunction, MyEqualFunction> prev;

    auto my_comp =
        [](const my_pair_t &e1, const my_pair_t &e2) { return e1.second > e2.second; };

    std::priority_queue<my_pair_t,
                        my_container_t,
                        decltype(my_comp)>
        queue(my_comp);

    queue.push(make_pair(current, 0));
    
    while (!queue.empty())
    {
        point ppp = currPosition.first;
        currPosition = queue.top();
        point curr = currPosition.first;

        queue.pop();
        if (curr == end)
        {
            prev[curr] = ppp;
            path.push_back(curr);
            break;
        }

        bool flag = false;

        point newPosition;
        size_t cost;

        //queue.empty();

        if (visited.find(curr) == visited.end())
        {
            visited.emplace(curr);
        }
        else
        {
            continue;
        }

        uint16_t localCost = 10000;
        point x ;
        if (m.can_go(0, curr.first, curr.second) && visited.find(curr + moveIn(0)) == visited.end())
        {
            newPosition = curr + moveIn(0);

            //vec[newPosition.first][newPosition.second]
            cost = currPosition.second + m.cost(curr, 0);
            if (cellCost[newPosition] > cost)
            {
                cellCost[newPosition] = cost;
                //queue.push(make_pair(newPosition, cost));
                //queue[newPosition] = cost;
                flag = true;
            }
        }
        if (m.can_go(1, curr.first, curr.second) && visited.find(curr + moveIn(1)) == visited.end())
        {
            newPosition = curr + moveIn(1);
            cost = currPosition.second + m.cost(curr, 1);
            if (cellCost[newPosition] > cost)
            {
                cellCost[newPosition] = cost;
                //queue.push(make_pair(newPosition, cost));
                flag = true;
            }
        }
        if (m.can_go(2, curr.first, curr.second) && visited.find(curr + moveIn(2)) == visited.end())
        {
            newPosition = curr + moveIn(2);
            cost = currPosition.second + m.cost(curr, 2);
            if (cellCost[newPosition] > cost)
            {
                cellCost[newPosition] = cost;
                //queue.push(make_pair(newPosition, cost));
                flag = true;
            }
        }
        if (m.can_go(3, curr.first, curr.second) && visited.find(curr + moveIn(3)) == visited.end())
        {
            newPosition = curr + moveIn(3);
            cost = currPosition.second + m.cost(curr, 3);
            if (cellCost[newPosition] > cost)
            {
                cellCost[newPosition] = cost;
                //queue.push(make_pair(newPosition, cost));
                flag = true;
            }
        }
        if (flag)
        {
            point x = queue.top().first;
            prev[x] =  currPosition.first;

        cout <<"stroing , prev of "<< prev[x].first <<", "<<prev[x].second <<" is " << 
        x.first <<", "<<x.second << "cost " << queue.top().second << "\n";

            path.push_back(x);
        }
        else
        {
            path.pop_back();
        }

        
    }

    path.clear();
    point currPt = end;
    while (prev.find(currPt) != prev.end() || currPt != make_pair(0,0))
    {
        point prevPT = prev[currPt];
        cout <<"prev of "<< currPt.first <<", "<<currPt.second <<" is " << 
        prevPT.first <<", "<<prevPT.second <<"\n";
        path.push_back(currPt);
        currPt = prevPT;
    }

    m.print_maze_with_path(cout, path, true, false);
    
    return path;
}

path solve_tour(Maze &m, int rows, int cols)
{
    return list<point>();
}
