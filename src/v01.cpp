#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <vector>

#define MPE_POLY2TRI_IMPLEMENTATION
#define MPE_POLY2TRI_USE_DOUBLE

#include <boost/multiprecision/gmp.hpp>
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

std::vector<SimplePolygon> triangulate(const Polygon &instance) {
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

        std::vector<SimplePolygon> triangles;
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

std::vector<SimplePolygon> solve(const Polygon &instance) {
    return triangulate(instance);
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
