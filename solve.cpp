#include "maze.h"
#include "path.h"
#include<queue>
#include<vector>
#include<list>
#include<tuple>
#include<utility>
#include<iostream>
#include<climits>
#include<sstream>
#include<stack>
#include<unordered_set>
using namespace std;

path solve_dfs(Maze& m, int rows, int cols);
path solve_bfs(Maze& m, int rows, int cols);
path solve_dijkstra(Maze& m, int rows, int cols);
path solve_tour(Maze& m, int rows, int cols);


class MyHashFunction { 
public: 
    // id is returned as hash function 
    size_t operator()(const point& p) const
    { 
        return (hash<int>()(p.first)) ^ (hash<int>()(p.second)); 
    } 
}; 

class MyEqualFunction { 
public: 
    bool operator()(const point& p1, const point& p2) const
    { 
        return ((p1.first == p2.first) && (p1.second == p2.second));
    } 
};

bool getUnvisitedNeighbors(point& neighbor,
                           unordered_set<point, MyHashFunction, MyEqualFunction>& visited, 
                           point currPosition,
                           Maze& m)
{
    for(uint16_t i = 0 ; i < 4 ; i++)
    {
        neighbor = currPosition + moveIn(i);
        if(m.can_go(i,currPosition.first, currPosition.second) && 
           visited.find(neighbor) == visited.end())
        {
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv)
{
    if(argc != 4)
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

    if(opt == "-dfs")
    {
        cout << "\nSolved dfs" << endl;
        path p = solve_dfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);
    }


    if(opt == "-bfs")
    {
        cout << "\nSolved bfs" << endl;
        path p = solve_bfs(m, rows, cols);
        m.print_maze_with_path(cout, p, false, false);
    }

    if(opt == "-dij")
    {
        cout << "\nSolved dijkstra" << endl;
        path p = solve_dijkstra(m, rows, cols);
        m.print_maze_with_path(cout, p, true, false);
    }

    if(opt == "-tour")
    {
        cout << "\nSolved all courners tour" << endl;
        path p = solve_tour(m, rows, cols);
        m.print_maze_with_path(cout, p, true, true);
    }
    if(opt == "-basic")
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
    if(opt == "-advanced")
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

path solve_dfs(Maze& m, int rows, int cols)
{

    stack<point> stacklist;
    point current = make_pair(0,0);
    stacklist.push(current);

    int row = m.rows();
    int col = m.columns();
    point currPosition;
    point end = make_pair(row-1,col-1);
    bool flag = false;

    unordered_set<point, MyHashFunction, MyEqualFunction> visited;
    list<point> path;

    while(true)
    {
        currPosition = stacklist.top();
        if(currPosition == end)
            break;
        point newPosition;

        if(visited.find(currPosition) == visited.end())
            visited.emplace(currPosition);

        bool found  = getUnvisitedNeighbors(newPosition,visited, currPosition,m);

        if(found)
        {
            stacklist.push(newPosition);
            path.push_back(newPosition);
        }
        else{
            stacklist.pop();
            path.pop_back();
        }
    }
    return path;
}

path solve_bfs(Maze& m, int rows, int cols)
{
    
    return list<point>();
}

path solve_dijkstra(Maze& m, int rows, int cols)
{
    return list<point>();
}

path solve_tour(Maze& m, int rows, int cols)
{
    return list<point>();
}
