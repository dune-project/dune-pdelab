// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
#ifndef DUNE_PDELAB_MULTITYPETREE_UTILITY_HH
#define DUNE_PDELAB_MULTITYPETREE_UTILITY_HH

#include <vector>
#include <iostream>
#include "multitypetree.hh"

namespace Dune {
  namespace PDELab {
    namespace MultiTypeTree {

      //!
      struct VisitingFunctor
      {
        void next() {}
        void down() {}
        void up() {}
        template<typename N> void enter_node(const N& n) {}
        template<typename N> void leave_node(const N& n) {}
        template<typename N> void visit_leaf(const N& n) {}
      };
      
      //!
      template<typename T, bool isleaf = T::isLeaf>
      struct MultiTypeTreeNodeVisitor;

      //!
      template<typename T, int n=T::CHILDREN, int i=0>
      struct MultiTypeTreeChildIterator;

#ifndef DOXYGEN
      // loop through childs
      template<typename T, int n, int i>
      struct MultiTypeTreeChildIterator
      {
        template<class F>
        static void visit_nodes (F & f, const T& t)
        {
          typedef typename T::template Child<i>::Type C;
          // visit child node
          MultiTypeTreeNodeVisitor<C>::visit_nodes(f,t.template getChild<i>());
          // go to next child
          f.next();
          MultiTypeTreeChildIterator<T,n,i+1>::visit_nodes(f,t);
        }
      };

      // behind last child (do nothing)... up again
      template<typename T, int n>
      struct MultiTypeTreeChildIterator<T,n,n>
      {
        template<class F>
        static void visit_nodes (F & f, const T& t) {}
      };

      // visit inner node
      template<typename T> 
      struct MultiTypeTreeNodeVisitor<T,false>
      {
        template<class F>
        static void visit_nodes (F & f, const T& t)
        {
          f.enter_node(t);
          f.down();
          MultiTypeTreeChildIterator<T>::visit_nodes(f,t);
          f.up();
          f.leave_node(t);
        }
      };
      
      // visit leaf node
      template<typename T> 
      struct MultiTypeTreeNodeVisitor<T,true>
      {
        template<class F>
        static void visit_nodes (F & f, const T& t)
        {
          f.visit_leaf(t);
        }
      };
#endif

      template<bool ignore_root=true, bool print_node=false, bool print_after=false>
      struct PrintPathFunctor : public VisitingFunctor {
        std::vector<int> path;
        std::ostream & os;
        PrintPathFunctor(std::ostream & os_) : os(os_) {
          if (!ignore_root)
            path.push_back(0);
        }
        void next() {
          path.back()++;
        }
        void down() {
          // and move down
          path.push_back(0);
        }
        void up() {
          // move back up
          path.pop_back();
        }
        template<typename P> void enter_node(const P& p) {
          // print
          if (print_node && !print_after)
          {
            os << "/";
            for (std::stack<std::string>::size_type i=0; i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
        }
        template<typename P> void leave_node(const P& p) {
          // and print
          if (print_node && print_after)
          {
            os << "/";
            for (std::stack<std::string>::size_type i=0; i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
        }
        template<typename P> void visit_leaf(const P& p) {
          os << "/";
          for (std::stack<std::string>::size_type i=0; i<path.size(); i++)
            os << path[i] << "/";
          os << std::endl;
        }
      };

      template<bool ignore_root=true, bool print_node=false, bool print_after=false>
      struct PrintPathFunctor2 : public VisitingFunctor {
        std::vector<int> path;
        std::ostream & os;
        PrintPathFunctor2(std::ostream & os_) : os(os_) {
            path.push_back(0);
        }
        template<typename N> void enter_node(const N& p) {
          // print
          if (print_node && !print_after)
          {
            os << "/";
            for (std::stack<std::string>::size_type i=(ignore_root?1:0);
                 i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
          // and move down
          path.push_back(0);
        }
        template<typename N> void leave_node(const N& p) {
          // move back up
          path.pop_back();
          // and print
          if (print_node && print_after)
          {
            os << "/";
            for (std::stack<std::string>::size_type i=(ignore_root?1:0);
                 i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
          // move to next node
          path.back()++;
        }
        template<typename N> void visit_leaf(const N& p) {
          os << "/";
          for (std::stack<std::string>::size_type i=(ignore_root?1:0);
               i<path.size(); i++)
            os << path[i] << "/";
          os << std::endl;
          // move to next node
          path.back()++;
        }
      };
      
      struct LeafCountFunctor : public VisitingFunctor {
        int i;
        LeafCountFunctor() : i(0) {}
        template<typename P> void visit_leaf(const P& p) {
          i++;
        }
      };

      struct NodeCountFunctor : public VisitingFunctor {
        int i;
        NodeCountFunctor() : i(0) {}
        template<typename P> void enter_node(const P& p) {
          i++;
        }
        template<typename P> void visit_leaf(const P& p) {
          i++;
        }
      };

      //!
      template<typename Functor, typename Tree>
      void ForEachNode(Functor& f, Tree& tree) {
          MultiTypeTreeNodeVisitor<Tree>::visit_nodes(f,tree);
      }

    } // end namespace MultiTypeTree

    /** \brief Count the number of nodes in a \ref MultiTypeTree
     *
     *  \tparam T A \ref MultiTypeTree
     *  \param  t An object of type T
     */
	template<typename T>
	int count_nodes (const T& t)
	{
      MultiTypeTree::NodeCountFunctor f;
      MultiTypeTree::ForEachNode(f,t);
	  return f.i;
	}

    /** \brief Count the number of leaves in a \ref MultiTypeTree
     *
     *  \tparam T A \ref MultiTypeTree
     *  \param  t An object of type T
     */
	template<typename T>
	int count_leaves (const T& t)
	{
      MultiTypeTree::LeafCountFunctor f;
      MultiTypeTree::ForEachNode(f,t);
	  return f.i;
	}

    /** \brief Print the paths leading to all subtrees in a \ref MultiTypeTree
     *
     *  \tparam T A \ref MultiTypeTree
     *  \param  t An object of type T
     */
    template<bool ignore_root = true,
             bool print_nodes = false,
             bool print_after = false,
             typename T>
	void print_paths (const T& t, std::ostream & s = std::cout)
	{
      MultiTypeTree::PrintPathFunctor<ignore_root, print_nodes, print_after> f(s);
      MultiTypeTree::ForEachNode(f,t);
	}

    template<typename T>
	void print_paths2 (const T& t, std::ostream & s = std::cout)
	{
      static const bool ignore_root = true;
      static const bool print_nodes = false;
      MultiTypeTree::PrintPathFunctor2<ignore_root, print_nodes> f(s);
      MultiTypeTree::ForEachNode(f,t);
	}

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_MULTITYPETREE_UTILITY_HH
