#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <list>
#include <optional>
#include <stdexcept>
#include <string>
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

SimplePolygon createRectangle(int xLeft, int yDown, int xRight, int yUp) {
    std::vector<Point> points;

    points.emplace_back(xLeft, yDown);
    points.emplace_back(xRight, yDown);
    points.emplace_back(xRight, yUp);
    points.emplace_back(xLeft, yUp);

    return {points.begin(), points.end()};
}

// These instances are small enough to be solved by hand
// Their solutions are implemented in code so that they can be verified using existing scripts
std::unordered_map<std::string, std::list<SimplePolygon>> solutionsByName{
        {"socg60",       {
                                 createRectangle(16, 592, 64, 624),
                                 createRectangle(16, 640, 64, 672),
                                 createRectangle(64, 560, 109, 704),
                                 createRectangle(112, 528, 144, 560),
                                 createRectangle(109, 592, 144, 672),
                                 createRectangle(112, 704, 144, 736),
                                 createRectangle(144, 480, 176, 800),
                                 createRectangle(176, 528, 192, 736),
                                 createRectangle(192, 480, 224, 800),
                                 createRectangle(224, 528, 256, 560),
                                 createRectangle(224, 592, 259, 672),
                                 createRectangle(224, 704, 256, 736),
                                 createRectangle(259, 560, 304, 704),
                                 createRectangle(304, 592, 352, 624),
                                 createRectangle(304, 640, 352, 672)
                         }},
        {"socg92",       {
                                 createRectangle(64, 176, 112, 208),
                                 createRectangle(64, 224, 112, 256),
                                 createRectangle(112, 144, 157, 288),
                                 createRectangle(160, 112, 192, 144),
                                 createRectangle(157, 176, 192, 256),
                                 createRectangle(160, 288, 192, 320),
                                 createRectangle(192, 64, 224, 384),
                                 createRectangle(224, 112, 240, 320),
                                 createRectangle(240, 64, 272, 384),
                                 createRectangle(272, 112, 304, 144),
                                 createRectangle(272, 288, 304, 320),
                                 createRectangle(272, 176, 330, 256),
                                 createRectangle(314, 112, 330, 144),
                                 createRectangle(314, 288, 330, 320),
                                 createRectangle(330, 64, 362, 384),
                                 createRectangle(362, 112, 378, 320),
                                 createRectangle(378, 64, 410, 384),
                                 createRectangle(410, 112, 442, 144),
                                 createRectangle(410, 176, 448, 256),
                                 createRectangle(410, 288, 442, 320),
                                 createRectangle(448, 144, 493, 288),
                                 createRectangle(493, 176, 538, 208),
                                 createRectangle(493, 224, 538, 256)
                         }},
        {"socg_fixed60", {
                                 createRectangle(240, 528, 288, 560),
                                 createRectangle(240, 576, 288, 608),
                                 createRectangle(288, 502, 320, 634),
                                 createRectangle(326, 464, 352, 496),
                                 createRectangle(320, 528, 352, 608),
                                 createRectangle(326, 640, 352, 672),
                                 createRectangle(352, 416, 384, 720),
                                 createRectangle(384, 464, 400, 672),
                                 createRectangle(400, 416, 432, 720),
                                 createRectangle(432, 464, 458, 496),
                                 createRectangle(432, 528, 464, 608),
                                 createRectangle(432, 640, 458, 672),
                                 createRectangle(464, 502, 496, 634),
                                 createRectangle(496, 528, 544, 560),
                                 createRectangle(496, 576, 544, 608)
                         }},
        {"socg_fixed92", {
                                 createRectangle(96, 256, 128, 320),
                                 createRectangle(96, 352, 128, 416),
                                 createRectangle(96, 544, 128, 608),
                                 createRectangle(96, 640, 128, 704),
                                 createRectangle(128, 224, 176, 448),
                                 createRectangle(128, 512, 176, 736),
                                 createRectangle(192, 160, 224, 208),
                                 createRectangle(176, 256, 224, 416),
                                 createRectangle(176, 544, 224, 704),
                                 createRectangle(192, 752, 224, 800),
                                 createRectangle(224, 128, 288, 832),
                                 createRectangle(288, 160, 320, 800),
                                 createRectangle(320, 128, 384, 832),
                                 createRectangle(384, 160, 416, 208),
                                 createRectangle(384, 256, 432, 416),
                                 createRectangle(384, 544, 432, 704),
                                 createRectangle(384, 752, 416, 800),
                                 createRectangle(432, 224, 480, 448),
                                 createRectangle(432, 512, 480, 736),
                                 createRectangle(480, 256, 512, 320),
                                 createRectangle(480, 352, 512, 416),
                                 createRectangle(480, 544, 512, 608),
                                 createRectangle(480, 640, 512, 704)
                         }}
};

// Taken from Boost
template<typename T>
inline void combineHashes(std::size_t &seed, const T &thing) {
    std::hash<T> hasher;
    seed ^= hasher(thing) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
    template<>
    struct hash<Point> {
        std::size_t operator()(const Point &point) const {
            std::size_t res = 0;

            combineHashes(res, CGAL::to_double(point.x()));
            combineHashes(res, CGAL::to_double(point.y()));

            return res;
        }
    };

    template<>
    struct hash<std::list<SimplePolygon>::iterator> {
        std::size_t operator()(const std::list<SimplePolygon>::iterator &it) const {
            return std::hash<SimplePolygon *>()(&*it);
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

std::optional<SimplePolygon> tryMergePolygons(const SimplePolygon &p1, const SimplePolygon &p2) {
    Polygon newPolygon;
    if (!CGAL::join(p1, p2, newPolygon)) {
        return std::nullopt;
    }

    if (newPolygon.has_holes()) {
        return std::nullopt;
    }

    const auto &boundary = newPolygon.outer_boundary();
    if (!boundary.is_convex()) {
        return std::nullopt;
    }

    return boundary;
}

std::optional<SimplePolygon> tryMergePolygons(const std::vector<SimplePolygon> &polygons) {
    std::vector<Polygon> newPolygons;
    CGAL::join(polygons.begin(), polygons.end(), std::back_inserter(newPolygons));

    if (newPolygons.size() != 1) {
        return std::nullopt;
    }

    const auto &polygon = newPolygons.at(0);
    if (polygon.has_holes()) {
        return std::nullopt;
    }

    const auto &boundary = polygon.outer_boundary();
    if (!boundary.is_convex()) {
        return std::nullopt;
    }

    return boundary;
}

void tryMergePolygonPairs(std::list<SimplePolygon> &polygons) {
    std::unordered_map<Point, std::unordered_set<std::list<SimplePolygon>::iterator>> polygonsByPoint;

    std::list<SimplePolygon>::iterator nextPolygon1;
    for (auto polygon1 = polygons.begin(); polygon1 != polygons.end(); polygon1 = nextPolygon1) {
        // We create a copy of the iterator because polygon1 may be invalid by the end of the iteration
        nextPolygon1 = polygon1;
        nextPolygon1++;

        for (const auto &point1 : *polygon1) {
            auto &set = polygonsByPoint[point1];
            if (set.empty()) {
                set.insert(polygon1);
                continue;
            }

            bool destroyedPolygon1 = false;

            for (auto polygon2 : set) {
                if (polygon1 == polygon2) {
                    continue;
                }

                auto newPolygon = tryMergePolygons(*polygon1, *polygon2);
                if (!newPolygon.has_value()) {
                    continue;
                }

                for (const auto &point : *polygon1) {
                    polygonsByPoint[point].erase(polygon1);
                }

                for (const auto &point : *polygon2) {
                    polygonsByPoint[point].erase(polygon2);
                }

                polygons.erase(polygon1);
                polygons.erase(polygon2);

                polygons.push_back(*newPolygon);

                destroyedPolygon1 = true;
                break;
            }

            if (destroyedPolygon1) {
                break;
            }

            set.insert(polygon1);
        }
    }
}

void tryMergePolygonsAll(std::list<SimplePolygon> &polygons) {
    std::unordered_map<Point, std::unordered_set<std::list<SimplePolygon>::iterator>> polygonsByPoint;

    std::list<SimplePolygon>::iterator nextPolygon1;
    for (auto polygon1 = polygons.begin(); polygon1 != polygons.end(); polygon1 = nextPolygon1) {
        // We create a copy of the iterator because polygon1 may be invalid by the end of the iteration
        nextPolygon1 = polygon1;
        nextPolygon1++;

        for (const auto &point1 : *polygon1) {
            auto &set = polygonsByPoint[point1];
            set.insert(polygon1);
            if (set.size() < 2) continue;

            std::vector<SimplePolygon> polys;
            for (auto poly : set) polys.push_back(*poly);
            auto newPolygon = tryMergePolygons(polys);
            if (!newPolygon.has_value()) {
                continue;
            }

            std::unordered_set<std::list<SimplePolygon>::iterator> copy(set.begin(), set.end());
            for (auto polygon : copy) {
                for (const auto &point : *polygon) {
                    polygonsByPoint[point].erase(polygon);
                }

                polygons.erase(polygon);
            }

            polygons.push_back(*newPolygon);
            break;
        }
    }
}

std::list<SimplePolygon> solve(const Polygon &instance, const std::string &name) {
    if (solutionsByName.contains(name)) {
        return solutionsByName[name];
    }

    std::list<SimplePolygon> polygons = triangulate(instance);
    tryMergePolygonsAll(polygons);
    tryMergePolygonPairs(polygons);

    return polygons;
}

SimplePolygon jsonToPolygon(const json &data) {
    std::vector<Point> points;

    for (const auto &p : data) {
        auto x = p["x"].get<std::int64_t>();
        auto y = p["y"].get<std::int64_t>();
        points.emplace_back(x, y);
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

    return {
            {"num", num},
            {"den", den}
    };
}

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    json input;
    if (argc < 2) {
        input = json::parse(std::cin);
    } else {
        std::ifstream fileStream(argv[1]);
        input = json::parse(fileStream);
    }

    SimplePolygon outerBoundary = jsonToPolygon(input["outer_boundary"]);

    std::vector<SimplePolygon> holes;
    for (const auto &hole : input["holes"]) {
        holes.emplace_back(jsonToPolygon(hole));
    }

    Polygon instance(outerBoundary, holes.begin(), holes.end());
    auto name = input["name"].get<std::string>();

    auto polygons = solve(instance, name);

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
            {"instance", name},
            {"polygons", outputPolygons}
    };

    std::cout << output << std::endl;
    std::cerr << polygons.size() << std::endl;

    auto end = std::chrono::high_resolution_clock::now();;
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cerr << duration.count() << "ms" << std::endl;

    return 0;
}
