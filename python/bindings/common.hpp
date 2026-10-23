#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace nb = nanobind;

using Value = int64_t;
using NextId = uint64_t;
using Coord = double;
using IdArray = nb::ndarray<nb::numpy, Value>;

inline void check_coords_shape(const nb::ndarray<nb::ro, Coord, nb::shape<-1, -1>>& coords, size_t dims, const char* label) {
    if (coords.ndim() != 2) 
        throw invalid_argument(string(label) + ": coords must be a 2D array");
    if (static_cast<size_t>(coords.shape(1)) != dims) 
        throw invalid_argument(string(label) + ": expected " + to_string(dims) + " columns in coords");
}

inline void check_corner_shape(const nb::ndarray<nb::ro, Coord, nb::shape<-1>>& corner, size_t dims, const char* label) {
    if (static_cast<size_t>(corner.shape(0)) != dims)
        throw invalid_argument(string(label) + ": expected length " + to_string(dims));
}

inline void check_ids_length(const nb::ndarray<nb::ro, Value, nb::shape<-1>>& ids, size_t n, const char* label) {
    if (static_cast<size_t>(ids.shape(0)) != n) 
        throw invalid_argument(string(label) + ": ids length must match number of rows in coords");
}

inline Value next_id_to_value(NextId next_id) {
    if (next_id > static_cast<NextId>(std::numeric_limits<Value>::max()))
        throw overflow_error("automatic id exhausted");
    return static_cast<Value>(next_id);
}

inline void advance_next_id_past(NextId& next_id, Value id) {
    if (id >= 0) {
        const NextId as_u = static_cast<NextId>(id);
        if (as_u >= next_id)
            next_id = as_u + 1;
    }
}

inline Value resolve_insert_id(const nb::object& id_obj, NextId& next_id) {
    if (id_obj.is_none()) {
        const Value result = next_id_to_value(next_id);
        ++next_id;
        return result;
    }
    const Value id = nb::cast<Value>(id_obj);
    advance_next_id_past(next_id, id);
    return id;
}

inline void update_next_id_after_bulk(const nb::object& ids_obj, const vector<Value>& id_values, NextId& next_id) {
    if (ids_obj.is_none()) {
        next_id = id_values.size();
        return;
    }
    NextId max_nonneg = 0;
    bool found_nonneg = false;
    for (Value v : id_values) {
        if (v >= 0) {
            const NextId as_u = static_cast<NextId>(v);
            if (!found_nonneg || as_u > max_nonneg) {
                max_nonneg = as_u;
                found_nonneg = true;
            }
        }
    }
    next_id = found_nonneg ? max_nonneg + 1 : 0;
}

inline vector<Value> resolve_ids(const nb::object& ids_obj, size_t n, const char* label) {
    vector<Value> ids(n);
    if (ids_obj.is_none()) {
        for (size_t i = 0; i < n; ++i)
            ids[i] = static_cast<Value>(i);
        return ids;
    }
    const auto arr = nb::cast<nb::ndarray<nb::ro, Value, nb::shape<-1>>>(ids_obj);
    check_ids_length(arr, n, label);
    for (size_t i = 0; i < n; ++i)
        ids[i] = arr(static_cast<int>(i));
    return ids;
}

inline IdArray empty_array() {
    return nb::ndarray<nb::numpy, Value>(nullptr, {0}, nb::handle());
}

inline IdArray to_array(vector<Value> values) {
    const size_t n = values.size();
    if (n == 0) 
        return empty_array();
    auto* data = new Value[n];
    copy(values.begin(), values.end(), data);
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<Value*>(p); });
    return nb::ndarray<nb::numpy, Value>(data, {n}, owner);
}

inline nb::dict stats_dict(size_t height, size_t count, size_t num_leaves, size_t num_internals, size_t size_bytes, size_t extra = 0, const char* extra_key = nullptr) {
    nb::dict d;
    d["height"] = height;
    d["count"] = count;
    d["num_leaves"] = num_leaves;
    d["num_internals"] = num_internals;
    d["size_bytes"] = size_bytes;
    if (extra_key != nullptr) 
        d[extra_key] = extra;
    return d;
}
