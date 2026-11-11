#pragma once

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct DIMACSInstance {
  size_t num_nodes;
  std::vector<std::pair<int, int>> interference_edges;
};

inline DIMACSInstance ParseDIMACS(const std::string &filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open DIMACS file: " + filepath);
  }

  DIMACSInstance instance{0, {}};
  std::string line;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == 'c') {
      continue; // Skip comments and empty lines
    }

    if (line[0] == 'p') {
      std::istringstream iss(line);
      std::string p_char, format;
      size_t num_edges;

      // Format: p edge <num_nodes> <num_edges>
      if (iss >> p_char >> format >> instance.num_nodes >> num_edges) {
        // Pre-allocate the exact size to prevent dynamic resizing during
        // parsing
        instance.interference_edges.reserve(num_edges);
      }
    } else if (line[0] == 'e') {
      std::istringstream iss(line);
      std::string e_char;
      int u, v;

      // Format: e <u_node> <v_node>
      if (iss >> e_char >> u >> v) {
        // DIMACS files are 1-indexed. Subtract 1 to match your 0-indexed
        // Grid1D.
        instance.interference_edges.push_back({u - 1, v - 1});
      }
    }
  }

  return instance;
}
