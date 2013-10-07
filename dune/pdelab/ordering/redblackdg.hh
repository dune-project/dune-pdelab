// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_REDBLACKDG_HH
#define DUNE_PDELAB_ORDERING_REDBLACKDG_HH

#include <vector>
#include <tuple>

#include <dune/common/exceptions.hh>

namespace Dune {
  namespace PDELab {

    template<typename GV, typename size_type = std::size_t>
    std::tuple<std::vector<size_type>,size_type> redBlackDGOrdering(const GV& gv)
    {
      enum class Color { unknown, red, black };
      std::vector<Color> color(gv.size(0),Color::unknown);
      std::vector<size_type> uncolored_neighbors;
      for (auto it = gv.template begin<0>(), end = gv.template end<0>(); it != end; ++it)
        {
          uncolored_neighbors.clear();
          auto& e = *it;
          size_type index = gv.indexSet().index(e);
          if (color[index] == Color::unknown)
            {
              // This cell hasn't been visited, look at its neighbors to derive its color
              Color neighbor_color = Color::unknown;
              for (auto iit = gv.ibegin(e), iend = gv.iend(e); iit != iend; ++iit)
                {
                  // skip boundary intersections
                  auto& is = *iit;
                  if (!is.neighbor())
                    continue;
                  auto oep = is.outside();
                  auto& oe = *oep;
                  size_type outside_index = gv.indexSet().index(oe);
                  Color c = color[outside_index];
                  if (c == Color::unknown)
                    uncolored_neighbors.push_back(outside_index);
                  if (neighbor_color == Color::unknown)
                    {
                      neighbor_color = c;
                    }
                  else if (c != Color::unknown)
                    {
                      if (c != neighbor_color)
                        DUNE_THROW(Exception,"Cannot color grid, inconsistent color in neighbors of cell " << index);
                    }

                }
              switch (neighbor_color)
                {
                case Color::unknown:
                  // we can pick a color
                  color[index] = Color::red;
                  neighbor_color = Color::black;
                  break;
                case Color::red:
                  color[index] = Color::black;
                  break;
                case Color::black:
                  color[index] = Color::red;
                  break;
                default:
                  DUNE_THROW(Exception,"This really shouldn't happen!");
                }

              // propagate color to neighbors
              for (auto outside_index : uncolored_neighbors)
                color[outside_index] = neighbor_color;
            }
          else
            {
              // Propagate / check color on neighbors
              Color neighbor_color = color[index] == Color::red ? Color::black : Color::red;

              for (auto iit = gv.ibegin(e), iend = gv.iend(e); iit != iend; ++iit)
                {
                  // skip boundary intersections
                  auto& is = *iit;
                  if (!is.neighbor())
                    continue;
                  auto oep = is.outside();
                  auto& oe = *oep;
                  size_type outside_index = gv.indexSet().index(oe);
                  Color c = color[outside_index];
                  if (c == Color::unknown)
                    color[outside_index] = neighbor_color;
                  else if (c != neighbor_color)
                    DUNE_THROW(Exception,"Cannot color grid, inconsistent color in neighbors of cell " << index);
                }
            }
        }

      std::vector<size_type> permutation(gv.size(0));

      size_type index = 0;
      for (size_type i = 0; i < permutation.size(); ++i)
        if (color[i] == Color::red)
          permutation[i] = index++;

      size_type black_offset = index;

      for (size_type i = 0; i < permutation.size(); ++i)
        if (color[i] == Color::black)
          permutation[i] = index++;

      return std::make_tuple(std::move(permutation),black_offset);

    }

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_REDBLACKDG_HH
