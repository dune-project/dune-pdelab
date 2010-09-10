// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
#ifndef DUNE_PDELAB_MULTITYPETREE_UTILITY_HH
#define DUNE_PDELAB_MULTITYPETREE_UTILITY_HH

#include <vector>
#include <iostream>
#include "multitypetree.hh"

namespace Dune {
  namespace PDELab {
    namespace MultiTypeTree {

      /** \addtogroup MultiTypeTreeTMP MultiTypeTree Template Meta Programs
       *  \ingroup MultiTypeTree
       *  \{
       */
      /**
         \interface VisitingFunctor
         \brief Interface of a visiting functor
         \see ForEachNode
       */
      struct VisitingFunctor
      {
        //! called before we move to a sibling node
        void next() {}
        //! called before we move a level down in the tree (to a set of childs)
        void down() {}
        //! called before we move a level up in the tree (back to the father node)
        void up() {}
        //! called when we have found a new inner node
        template<typename N> void enter_node(const N& n) {}
        //! called when we leave an inner node
        template<typename N> void leave_node(const N& n) {}
        //! called when we have found a new leaf node
        template<typename N> void visit_leaf(const N& n) {}
      };
      
      //! TMP implementing the "if" of the algortihm of ForEachNode
      //! \memberof ForEachNode
      template<typename T, bool isleaf = T::isLeaf>
      struct MultiTypeTreeNodeVisitor;

      //! TMP implementing the "for"-loop of the algortihm of ForEachNode
      //! \memberof ForEachNode
      template<typename T, int n=T::CHILDREN, int i=0>
      struct MultiTypeTreeChildIterator;

      // \ingroup MultiTypeTree
      /** \brief Loop through all nodes of a MultiTypeTree using depth-first and call a functor

          \param f     functor fulfilling the interface of VisitingFunctor
          \param tree  tree to loop

          The user has to prove a functor, fulfilling the interface of VisitingFunctor.
          This has functor implements a set, or subset, of hooks which allow to implement
          arbitrary operations n the tree.
          
          The function is implemented using as a Template Meta Program.

          The following psuedo code shows the alorithm an where the different hooks are applied:
          \code
          void ForEachNode(VisitingFunctor & f, const Node & tree)
          {
              if (n.is_leaf())
              {
                  // visit leaf
                  f.visit_leaf(n);
              }
              else
              {
                  // enter node
                  f.enter_node(n);
                  // down
                  f.down();
                  for (unsigned int i=0; i<n.size(); i++)
                  {
                      ForEachNode(f, n.get(i));
                      // next
                      f.next();
                  }
                  // up
                  f.up();
                  // leave node
                  f.leave_node(n);
              }
          }
          \endcode
      */
      template<typename Functor, typename Tree>
      void ForEachNode(Functor& f, Tree& tree) {
          MultiTypeTreeNodeVisitor<Tree>::visit_nodes(f,tree);
      }

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

      /** \brief VisitingFunctor to read all nodes in a MultiTypeTree
          \tparam store_root   should the root node be stored

          The functor notes all nodes an leafs and stores their path in an
          \code
          std::vector<int> path
          \endcode

          You can derive your own class and obtain the current path from the path member.
      */
      template<bool store_root=false>
      struct ReadPathFunctor : public VisitingFunctor {
        //! vector storing the current path
        std::vector<int> path;
        /** constructor */
        ReadPathFunctor() {
          if (store_root)
            path.push_back(0);
        }
        void next() {
          path.back()++;
        }
        void down() {
          path.push_back(0);
        }
        void up() {
          path.pop_back();
        }
      };

      /** \brief implement VisitingFunctor for print_paths

          \tparam print_root   should the root node be printed?
          \tparam print_node   should inner nodes be printed?
          \tparam print_after  should inner nodes be printed when entering the node (true), or when leaving the node (false)?
       */
      template<bool print_root=false, bool print_node=false, bool print_after=false>
      struct PrintPathFunctor : public ReadPathFunctor<print_root> {
        using ReadPathFunctor<print_root>::path;
        std::ostream & os;
        /** constructor */
        PrintPathFunctor(std::ostream & os_) : os(os_) {
        }
        template<typename P> void enter_node(const P& p) {
          // print
          if (print_node && !print_after)
          {
            os << "/";
            for (std::vector<int>::size_type i=0; i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
        }
        template<typename P> void leave_node(const P& p) {
          // and print
          if (print_node && print_after)
          {
            os << "/";
            for (std::vector<int>::size_type i=0; i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
        }
        template<typename P> void visit_leaf(const P& p) {
          os << "/";
          for (std::vector<int>::size_type i=0; i<path.size(); i++)
            os << path[i] << "/";
          os << std::endl;
        }
      };

#ifndef DOXYGEN
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
            for (std::vector<int>::size_type i=(ignore_root?1:0);
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
            for (std::vector<int>::size_type i=(ignore_root?1:0);
                 i<path.size(); i++)
              os << path[i] << "/";
            os << std::endl;
          }
          // move to next node
          path.back()++;
        }
        template<typename N> void visit_leaf(const N& p) {
          os << "/";
          for (std::vector<int>::size_type i=(ignore_root?1:0);
               i<path.size(); i++)
            os << path[i] << "/";
          os << std::endl;
          // move to next node
          path.back()++;
        }
      };
#endif

      //! implement VisitingFunctor for count_leaves
      struct LeafCountFunctor : public VisitingFunctor {
        int i;
        LeafCountFunctor() : i(0) {}
        template<typename P> void visit_leaf(const P& p) {
          i++;
        }
      };

      //! implement VisitingFunctor for count_nodes
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

      //! \}
    
    } // end namespace MultiTypeTree

    //! \ingroup MultiTypeTree
    //! \{
    
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
     *  \param  s output stream
     */
    template<bool print_root = false,
             bool print_nodes = false,
             bool print_after = false,
             typename T>
	void print_paths (const T& t, std::ostream & s = std::cout)
	{
      MultiTypeTree::PrintPathFunctor<print_root, print_nodes, print_after> f(s);
      MultiTypeTree::ForEachNode(f,t);
	}

#ifndef DOXYGEN
    template<typename T>
	void print_paths2 (const T& t, std::ostream & s = std::cout)
	{
      static const bool ignore_root = true;
      static const bool print_nodes = false;
      MultiTypeTree::PrintPathFunctor2<ignore_root, print_nodes> f(s);
      MultiTypeTree::ForEachNode(f,t);
	}
#endif

    //! \}
    
  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_MULTITYPETREE_UTILITY_HH
