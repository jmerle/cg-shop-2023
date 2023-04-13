#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <list>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#define MPE_POLY2TRI_IMPLEMENTATION
#define MPE_POLY2TRI_USE_DOUBLE

#include <boost/multiprecision/gmp.hpp>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <MPE_fastpoly2tri.h>
#include <nlohmann/json.hpp>

using Kernel = CGAL::Epeck;
using Point = CGAL::Point_2<Kernel>;
using Polygon = CGAL::Polygon_with_holes_2<Kernel>;
using SimplePolygon = CGAL::Polygon_2<Kernel>;

using json = nlohmann::json;

struct EdgeInfo {
    Point point1;
    Point point2;

    std::list<SimplePolygon>::iterator polygon;
    int idxPolygon;

    int idxPoint1;
    int idxPoint2;

    EdgeInfo(std::list<SimplePolygon>::iterator polygon, int idxPolygon, int idx1)
            : polygon(polygon), idxPolygon(idxPolygon) {
        const SimplePolygon poly = *polygon;

        int idx2 = idx1 < poly.vertices().size() - 1 ? idx1 + 1 : 0;

        const auto &p1 = poly.vertex(idx1);
        const auto &p2 = poly.vertex(idx2);

        if (p1 <= p2) {
            point1 = p1;
            point2 = p2;

            idxPoint1 = idx1;
            idxPoint2 = idx2;
        } else {
            point1 = p2;
            point2 = p1;

            idxPoint1 = idx2;
            idxPoint2 = idx1;
        }
    }

    bool operator==(const EdgeInfo &rhs) const {
        return point1 == rhs.point1 && point2 == rhs.point2;
    }

    bool operator!=(const EdgeInfo &rhs) const {
        return !(rhs == *this);
    }
};

// Taken from Boost
template<typename T>
inline void combineHashes(std::size_t &seed, const T &thing) {
    std::hash<T> hasher;
    seed ^= hasher(thing) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
    template<>
    struct hash<EdgeInfo> {
        std::size_t operator()(const EdgeInfo &edge) const {
            std::size_t res = 0;

            combineHashes(res, CGAL::to_double(edge.point1.x()));
            combineHashes(res, CGAL::to_double(edge.point1.y()));
            combineHashes(res, CGAL::to_double(edge.point2.x()));
            combineHashes(res, CGAL::to_double(edge.point2.y()));

            return res;
        }
    };
}

std::list<SimplePolygon> triangulate(const Polygon &instance) {
    std::uint32_t maxPointCount = instance.outer_boundary().size();
    for (const auto &hole : instance.holes()) {
        maxPointCount += hole.size();
    }

    std::size_t memoryRequired = MPE_PolyMemoryRequired(maxPointCount);
    void *memory = std::calloc(memoryRequired, 1);

    MPEPolyContext polyContext;
    if (MPE_PolyInitContext(&polyContext, memory, maxPointCount)) {
        for (const auto &p : instance.outer_boundary()) {
            MPEPolyPoint *point = MPE_PolyPushPoint(&polyContext);
            point->X = CGAL::to_double(p.x());
            point->Y = CGAL::to_double(p.y());
        }

        MPE_PolyAddEdge(&polyContext);

        for (const auto &hole : instance.holes()) {
            MPEPolyPoint *points = MPE_PolyPushPointArray(&polyContext, hole.size());

            for (int i = 0, iMax = hole.size(); i < iMax; i++) {
                const auto &p = hole[i];
                points[i].X = CGAL::to_double(p.x());
                points[i].Y = CGAL::to_double(p.y());
            }

            MPE_PolyAddHole(&polyContext);
        }

        MPE_PolyTriangulate(&polyContext);

        std::list<SimplePolygon> triangles;
        for (int i = 0; i < polyContext.TriangleCount; i++) {
            MPEPolyTriangle *triangle = polyContext.Triangles[i];

            std::vector<Point> points;
            for (const auto *p : triangle->Points) {
                points.emplace_back(p->X, p->Y);
            }

            triangles.emplace_back(points.begin(), points.end());
        }

        std::free(memory);
        return triangles;
    } else {
        std::free(memory);
        throw std::runtime_error("Could not initialize poly context");
    }
}

bool tryRemoveEdge(std::list<SimplePolygon> &polygons, const EdgeInfo &edge1, const EdgeInfo &edge2) {
    const auto &polygon1 = *edge1.polygon;
    const auto &polygon2 = *edge2.polygon;

    int size1 = polygon1.size();
    int size2 = polygon2.size();

    std::vector<Point> newPoints;
    newPoints.reserve(size1 + size2 - 2);

    int modifier1 = (edge1.idxPoint1 + 1) % size1 != edge1.idxPoint2 ? 1 : -1;
    for (int i = edge1.idxPoint1, j = 0;
         j < size1 - 1;
         i = (i + modifier1) >= 0 ? (i + modifier1) % size1 : size1 - 1, j++) {
        newPoints.push_back(polygon1.vertex(i));
    }

    int modifier2 = (edge2.idxPoint2 + 1) % size2 != edge2.idxPoint1 ? 1 : -1;
    for (int i = edge2.idxPoint2, j = 0;
         j < size2 - 1;
         i = (i + modifier2) >= 0 ? (i + modifier2) % size2 : size2 - 1, j++) {
        newPoints.push_back(polygon2.vertex(i));
    }

    SimplePolygon newPolygon(newPoints.begin(), newPoints.end());
    if (!newPolygon.is_convex()) {
        return false;
    }

    // The code above gets the new polygon correct, but with the wrong point order
    // The is_convex() check still works, but the area is sometimes negative, which the verifier doesn't like
    // CGAL::join() gets the point order right, but it's a lot slower, so the code above serve as a fast filter
    Polygon newPolygonSorted;
    CGAL::join(polygon1, polygon2, newPolygonSorted);

    polygons.erase(edge1.polygon);
    polygons.erase(edge2.polygon);

    polygons.push_back(newPolygonSorted.outer_boundary());
    return true;
}

std::list<SimplePolygon> solve(const Polygon &instance) {
    std::list<SimplePolygon> polygons = triangulate(instance);

    std::unordered_map<EdgeInfo, EdgeInfo> edgeMap;

    int i = 0;
    std::list<SimplePolygon>::iterator nextPolygon;
    for (auto polygon = polygons.begin(); polygon != polygons.end(); polygon = nextPolygon, i++) {
        nextPolygon = polygon;
        nextPolygon++;

        for (int j = 0, jMax = (*polygon).vertices().size(); j < jMax; j++) {
            EdgeInfo edge1(polygon, i, j);

            auto item = edgeMap.find(edge1);
            if (item == edgeMap.end()) {
                edgeMap.insert({edge1, edge1});
                continue;
            }

            const auto &edge2 = item->second;
            if (tryRemoveEdge(polygons, edge1, edge2)) {
                for (auto it = edgeMap.begin(); it != edgeMap.end();) {
                    if (it->second.idxPolygon == i || it->second.idxPolygon == edge2.idxPolygon) {
                        it = edgeMap.erase(it);
                    } else {
                        it++;
                    }
                }
                break;
            }
        }
    }

    return polygons;
}

Kernel::FT intToCGAL(std::int64_t value) {
    double low32 = value & 0xffffffff;
    double high32 = double(value >> 32) * 4294967296.0;
    return Kernel::FT(low32) + Kernel::FT(high32);
}

SimplePolygon jsonToPolygon(const json &data) {
    std::vector<Point> points;

    for (const auto &p : data) {
        auto x = p["x"].get<std::int64_t>();
        auto y = p["y"].get<std::int64_t>();
        points.emplace_back(intToCGAL(x), intToCGAL(y));
    }

    return {points.begin(), points.end()};
}

json coordinateToJSON(const Kernel::FT &coord) {
    const auto &value = coord.exact();

    auto num = numerator(value).convert_to<long long>();
    auto den = denominator(value).convert_to<long long>();

    if (den == 1) {
        return num;
    }

    return json{
            {"num", num},
            {"den", den}
    };
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    json input = json::parse(std::cin);

    SimplePolygon outerBoundary = jsonToPolygon(input["outer_boundary"]);

    std::vector<SimplePolygon> holes;
    for (const auto &hole : input["holes"]) {
        holes.emplace_back(jsonToPolygon(hole));
    }

    Polygon instance(outerBoundary, holes.begin(), holes.end());

    auto polygons = solve(instance);

    json outputPolygons = json::array();

    for (const auto &polygon : polygons) {
        json points = json::array();

        for (const auto &p : polygon.container()) {
            points.emplace_back(json{
                    {"x", coordinateToJSON(p.x())},
                    {"y", coordinateToJSON(p.y())}
            });
        }

        outputPolygons.emplace_back(points);
    }

    json output = {
            {"type",     "CGSHOP2023_Solution"},
            {"instance", input["name"]},
            {"polygons", outputPolygons}
    };

    std::cerr << polygons.size() << std::endl;
    std::cout << output << std::endl;

    return 0;
}
