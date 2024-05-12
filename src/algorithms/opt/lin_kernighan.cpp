#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <random>

#define NUMBER long long


#define EXTERN_DLL_EXPORT __attribute__((visibility("default")))
template<typename T>
class PointerTo1D {
private:
    T *data;
    const unsigned long _size;
public:
    PointerTo1D(T *data, NUMBER size) : data(data), _size(size) {}

    T operator[](unsigned long i) {
        return data[i];
    }
};

template<typename T>
class Array2DFrom1D {
private:
    std::vector<T> data;
    const unsigned long _size;

public:
    unsigned long size() const {
        return _size;
    }

    explicit Array2DFrom1D(NUMBER size) : _size(size) {
        data.resize(size * size);
    }

    Array2DFrom1D(std::vector<T> data, NUMBER size) : data(data), _size(size) {
        if (data.size() != size * size) {
            throw std::invalid_argument("Data size does not match size*size");
        }
    }
//
//    NUMBER &operator[](NUMBER i, NUMBER j) {
//        return data[i * size + j];
//    }

    PointerTo1D<NUMBER> operator[](NUMBER i) const {

        return PointerTo1D<NUMBER>((NUMBER *) &data[i * _size], _size);
    }

};


inline std::vector<NUMBER> nearest_neighbour_solution(const Array2DFrom1D<NUMBER> &dist_matrix) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, (int) (dist_matrix.size()) - 1);

    NUMBER node = dis(gen);
    std::vector<NUMBER> result = {node};

    std::vector<NUMBER> nodes_to_visit(dist_matrix.size());
    std::iota(nodes_to_visit.begin(), nodes_to_visit.end(), 0); // Fill with 0, 1, ..., size-1
    nodes_to_visit.erase(std::remove(nodes_to_visit.begin(), nodes_to_visit.end(), node),
                         nodes_to_visit.end()); // Remove initial node

    while (!nodes_to_visit.empty()) {
        auto nearest_node = std::min_element(nodes_to_visit.begin(), nodes_to_visit.end(),
                                             [&](NUMBER i, NUMBER j) {
                                                 return dist_matrix[node][i] < dist_matrix[node][j];
                                             });
        node = *nearest_node;
        nodes_to_visit.erase(std::remove(nodes_to_visit.begin(), nodes_to_visit.end(), node), nodes_to_visit.end());
        result.push_back(node);
    }

    return result;
}


struct Element {
public:
    NUMBER node;
    NUMBER cost;

    Element(NUMBER node, NUMBER cost) : node(node), cost(cost) {}

};

class Path {
public:
    std::vector<NUMBER> path;
    NUMBER cost;

    explicit Path(std::vector<NUMBER> path, NUMBER cost = 0) : path((path)), cost(cost) {}

    bool operator==(const Path &other) const {
        if (other.path.size() != path.size()) {
            return false;
        }
        for (NUMBER i = 0; i < path.size(); ++i) {
            if (path[i] != other.path[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator<(const Path &other) const {
        // Define what it means for one Path to be less than another
        return cost < other.cost;
    }
};

class Tour {
public:
    std::vector<NUMBER> tour_array;
    std::unordered_map<NUMBER, NUMBER> tour_indexes;
    NUMBER length;
    std::set<std::pair<NUMBER, NUMBER>> edges;

    Tour(std::vector<NUMBER> tour, Array2DFrom1D<NUMBER> data) : tour_array(tour) {
        for (NUMBER i = 0; i < tour.size(); ++i) {
            tour_indexes[tour[i]] = i;
        }
        length = cost(data);
        for (NUMBER i = 0; i < size(); ++i) {
            NUMBER a = tour_array[i - 1 >= 0 ? i - 1 : size() - 1];
            NUMBER b = tour_array[i];
            if (a < b) {
                edges.insert(std::make_pair(a, b));
            } else {
                edges.insert(std::make_pair(b, a));
            }
        }
    }

    NUMBER size() const {
        return tour_array.size();
    }

    NUMBER cost(const Array2DFrom1D<NUMBER> &data) {
        NUMBER cost = 0;
        for (NUMBER i = 0; i < size(); ++i) {
            cost += data[tour_array[i - 1 >= 0 ? i - 1 : tour_array.size() - 1]][tour_array[i]];
        }
        return cost;
    }

    NUMBER operator[](NUMBER index) const {
        return tour_array[index];
    }

    bool contains(const std::pair<NUMBER, NUMBER> &edge) const {
        return edges.find(edge) != edges.end() || edges.find(std::make_pair(edge.second, edge.first)) != edges.end();
    }

    NUMBER index(NUMBER node) const {
        return tour_indexes.at(node);
    }

    std::pair<NUMBER, NUMBER> around(NUMBER node) const {
        NUMBER index = this->index(node);
        NUMBER prev_child = index - 1 >= 0 ? index - 1 : size() - 1;
        NUMBER next_child = index + 1 < size() ? index + 1 : 0;
        return std::make_pair((*this)[prev_child], (*this)[next_child]);
    }

    std::pair<bool, std::vector<NUMBER>>
    generate(std::set<std::pair<NUMBER, NUMBER>> broken, std::set<std::pair<NUMBER, NUMBER>> joined) {
        std::set<std::pair<NUMBER, NUMBER>> edges = this->edges;
        for (auto &b: broken) {
            edges.erase(b);
        }
        for (auto &j: joined) {
            edges.insert(j);
        }

        if (edges.size() < this->size()) {
//            std::cerr << "Not enough edges to form a tour." << std::endl;
            return {false, {}};
        }
        std::unordered_map<NUMBER, NUMBER> successors = _build_successors(edges);
        if (successors.size() < this->size()) {
//            std::cerr << "Not enough successors to form a tour." << std::endl;
            return {false, {}};
        }

        NUMBER succ = successors[0];
        std::vector<NUMBER> new_tour = {0};
        std::set<NUMBER> visited = {0};

        while (visited.find(succ) == visited.end()) {
            visited.insert(succ);
            new_tour.push_back(succ);
            succ = successors[succ];
        }

        return {new_tour.size() == this->size(), new_tour};
    }

private:
    std::unordered_map<NUMBER, NUMBER> _build_successors(std::set<std::pair<NUMBER, NUMBER>> edges) {
        std::unordered_map<NUMBER, NUMBER> successors = {};
        NUMBER node = 0;


        while (!edges.empty()) {
            NUMBER first, second;
            auto it = edges.begin();
            while (it != edges.end()) {
                first = it->first;
                second = it->second;
                if (first == node) {
                    successors[node] = second;
                    node = second;
                    break;
                } else if (second == node) {
                    successors[node] = first;
                    node = first;
                    break;
                } else {
                    ++it;
                }
            }
            edges.erase(std::make_pair(first, second));
        }

        return successors;
    }
};


class LinKernighan {
public:
    const Array2DFrom1D<NUMBER> data;
    NUMBER startPoint;
    std::vector<NUMBER> heuristicPath;
    std::set<Path> paths;
    NUMBER cost;
    std::map<NUMBER, std::vector<NUMBER>> neighbours;

    LinKernighan(const Array2DFrom1D<NUMBER> &data, NUMBER startPoint,
                 const std::vector<NUMBER> &startPath = {}) :
            data(data), startPoint(startPoint), heuristicPath(startPath) {

        if (heuristicPath.empty()) {
            // Initialize heuristicPath with nearest neighbour solution
            heuristicPath = nearest_neighbour_solution(data);
        }
        cost = 0;
        for (auto i: heuristicPath) {
            cost += data[i][(i + 1) % heuristicPath.size()];
        }
        initializeNeighbours();
    }

    void initializeNeighbours() {
        for (auto i: heuristicPath) {
            std::vector<NUMBER> neighbours_;
            for (NUMBER j = 0; j < data.size(); ++j) {
                if (i != j) {
                    neighbours_.push_back(j);
                }
            }
            this->neighbours[i] = neighbours_;

        }
    }

    std::vector<Element> solve() {
        bool improve = true;
        while (improve) {
            improve = this->improve();
            auto g = Path(heuristicPath, cost);
            paths.insert(g);
        }

        auto bestPath = std::min_element(paths.begin(), paths.end());
        std::vector<Element> result;
        for (size_t i = 0; i < bestPath->path.size(); ++i) {
            NUMBER current = bestPath->path[i];
            NUMBER next = bestPath->path[(i + 1) % bestPath->path.size()];
            NUMBER cost = data[current][next];
            result.emplace_back(current, cost);
        }
        return result;
    }

    std::vector<std::pair<NUMBER, std::pair<NUMBER, NUMBER>>>
    closest(NUMBER t2i, Tour &tour, NUMBER gain, std::set<std::pair<NUMBER, NUMBER>> &broken,
            std::set<std::pair<NUMBER, NUMBER>> &joined) {
        std::map<NUMBER, std::pair<NUMBER, NUMBER>> neighboursMap;

        // Create the neighbours of t_2i
        for (NUMBER node: neighbours[t2i]) {
            std::pair<NUMBER, NUMBER> yi = (t2i < node) ? std::make_pair(t2i, node) : std::make_pair(node, t2i);
            NUMBER Gi = gain - data[t2i][node];

            // Continue if the conditions are not met
            if (Gi <= 0 || broken.find(yi) != broken.end() || tour.contains(yi)) {
                continue;
            }
            auto aroundIT = tour.around(node);
            auto aroundVector = std::vector<NUMBER>{aroundIT.first, aroundIT.second};
            for (NUMBER succ: aroundVector) {
                std::pair<NUMBER, NUMBER> xi = (node < succ) ? std::make_pair(node, succ) : std::make_pair(succ, node);
                if (broken.find(xi) != broken.end() || joined.find(xi) != joined.end()) {
                    continue;
                }

                NUMBER diff = data[node][succ] - data[t2i][node];
                auto it = neighboursMap.find(node);
                if (it != neighboursMap.end() && diff > it->second.first) {
                    it->second = std::make_pair(diff, it->second.second);
                } else {
                    neighboursMap[node] = std::make_pair(diff, Gi);
                }
            }
        }

        // Sort neighbours by potential improvement
        std::vector<std::pair<NUMBER, std::pair<NUMBER, NUMBER>>> sortedNeighbours;
        for (const auto &[node, values]: neighboursMap) {
            sortedNeighbours.emplace_back(node, std::make_pair(values.first, values.second));
        }

        std::sort(sortedNeighbours.begin(), sortedNeighbours.end(), [](const auto &a, const auto &b) {
            auto g = a.second.first;
            auto g2 = b.second.first;
            return g > g2;
        });

        return sortedNeighbours;
    }

    bool improve() {
        Tour tour(heuristicPath, data);
        for (NUMBER t1: heuristicPath) {
            NUMBER before, after;
            std::tie(before, after) = tour.around(t1);
            for (NUMBER t2: {before, after}) {
                std::set<std::pair<NUMBER, NUMBER>> broken = {
                        (t1 < t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1)};
                NUMBER gain = data[t1][t2];

                auto set_ = std::set<std::pair<NUMBER, NUMBER>>{};
                auto close = closest(t2, tour, gain, broken, set_);

                NUMBER tries = 5;

                for (const auto &[t3, gi_pair]: close) {
                    NUMBER Gi = gi_pair.second;

                    if (t3 == before || t3 == after) {
                        continue;
                    }

                    std::set<std::pair<NUMBER, NUMBER>> joined = {
                            (t2 < t3) ? std::make_pair(t2, t3) : std::make_pair(t3, t2)};

                    if (choice_remove(tour, t1, t3, Gi, broken, joined)) {
                        return true;
                    }

                    tries--;
                    if (tries <= 0) break;
                }
            }
        }
        return false;
    }

    bool choice_remove(Tour &tour, NUMBER t1, NUMBER last, NUMBER gain, std::set<std::pair<NUMBER, NUMBER>> &broken,
                       std::set<std::pair<NUMBER, NUMBER>> &joined) {
        std::vector<NUMBER> around;
        if (broken.size() == 4) {
            NUMBER pred, succ;
            std::tie(pred, succ) = tour.around(last);
            if (this->data[pred][last] > this->data[succ][last]) {
                around = {pred};
            } else {
                around = {succ};
            }
        } else {
            auto a = tour.around(last);
            around = {a.first, a.second};
        }
        for (NUMBER t2i: around) {
            std::pair<NUMBER, NUMBER> xi = (last < t2i) ? std::make_pair(last, t2i) : std::make_pair(t2i, last);
            NUMBER Gi = gain + data[last][t2i];

            if (joined.find(xi) != joined.end() || broken.find(xi) != broken.end()) {
                continue;
            }

            std::set<std::pair<NUMBER, NUMBER>> added(joined);
            std::set<std::pair<NUMBER, NUMBER>> removed(broken);

            removed.insert(xi);
            added.insert((t2i < t1) ? std::make_pair(t2i, t1) : std::make_pair(t1, t2i));
            NUMBER relink = Gi - data[t2i][t1];

            auto [isTour, newTour] = tour.generate(removed, added);
            if (!isTour && added.size() > 2) {
                continue;
            }

            auto existsPath = paths.find(Path(newTour, cost));
            if (existsPath != paths.end() && existsPath->operator==(Path(newTour, cost))) {
                return false;
            }

            if (isTour && relink > 0) {
                heuristicPath = newTour;
                cost -= relink;
                return true;
            } else {
                auto ans = choice_add(tour, t1, t2i, Gi, removed, joined);
                if (broken.size() == 2 && ans) {
                    return true;
                }
                return ans;
            }
        }

        return false;
    }

    bool choice_add(Tour &tour, NUMBER t1, NUMBER t2i, NUMBER gain, std::set<std::pair<NUMBER, NUMBER>> &broken,
                    std::set<std::pair<NUMBER, NUMBER>> &joined) {
        auto ordered = closest(t2i, tour, gain, broken, joined);

        NUMBER top = (broken.size() == 2) ? 5 : 1;

        for (const auto &[node, gi_pair]: ordered) {
            NUMBER Gi = std::get<1>(gi_pair);
            std::pair<NUMBER, NUMBER> yi = (t2i < node) ? std::make_pair(t2i, node) : std::make_pair(node, t2i);
            std::set<std::pair<NUMBER, NUMBER>> newJoined(joined);
            newJoined.insert(yi);

            if (choice_remove(tour, t1, node, Gi, broken, newJoined)) {
                return true;
            }

            if (--top == 0) {
                return false;
            }
        }

        return false;
    }

};

extern "C" {
EXTERN_DLL_EXPORT int solve(int size, NUMBER *data, NUMBER *start_data, NUMBER *result) {
    /**
     * size - size of data (horizontal and vertical)
     * data - 1D array of size*size
     * start_data - 1D array of size (start path), can be NULL,
     * result - 1D array of size
     */
    try {
        std::vector<NUMBER> test;
        test.resize(size * size);
        for (int i = 0; i < size * size; ++i) {
            test[i] = data[i];
        }
        auto data_ = Array2DFrom1D<NUMBER>(test, size);
        if (start_data != nullptr && start_data[0] != -1) {
            std::vector<NUMBER> start;
            start.resize(size);
            for (int i = 0; i < size; i++) {
                start[i] = start_data[i];
            }
            LinKernighan lk = LinKernighan(data_, 0, start);
            auto ans = lk.solve();
            for (int i = 0; i < size; i++) {
                result[i] = ans[i].node;
            }
        } else {
            LinKernighan lk = LinKernighan(data_, 0);
            auto ans = lk.solve();
            for (int i = 0; i < size; i++) {
                result[i] = ans[i].node;
            }
        }

        return 0;
    }
    catch (std::exception e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
}
