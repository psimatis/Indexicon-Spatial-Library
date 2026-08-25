#include <nanobind/nanobind.h>

namespace nb = nanobind;

void bind_rtree_point(nb::module_& m);
void bind_rtree_mbr(nb::module_& m);
void bind_quadtree(nb::module_& m);
void bind_mxcif_quadtree(nb::module_& m);
void bind_octtree(nb::module_& m);
void bind_kdtree(nb::module_& m);

NB_MODULE(_core, m) {
    m.doc() = "Indexicon spatial index bindings";
    m.attr("__indexicon_cpp_version__") = INDEXICON_CPP_VERSION;

    bind_rtree_point(m);
    bind_rtree_mbr(m);
    bind_quadtree(m);
    bind_mxcif_quadtree(m);
    bind_octtree(m);
    bind_kdtree(m);
}
