#ifndef DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH
#define DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

class GnuplotGraph
{
public:
  GnuplotGraph(const std::string &filename_)
    : mode(command)
    , filename(filename_)
    , stream(filename.c_str())
  {}

  void addCommand(const std::string &cmd)
  {
    commandMode();
    stream << cmd << std::endl;
  }

  void addPlot(const std::string &plotstuff)
  {
    plotMode();
    stream << plotDelim << plotstuff;
    plotDelim = ", \\\n     ";
  }

  ~GnuplotGraph()
  {
    commandMode();
    stream.close();
    std::ostringstream command;
    command << "gnuplot " << filename;
    std::system(command.str().c_str());
  }

private:
  void commandMode()
  {
    switch(mode) {
    case command:
      break;
    case plot:
      stream << std::endl;
      break;
    }
    mode = command;
  }
  
  void plotMode()
  {
    switch(mode) {
    case command:
      plotDelim = "plot \\\n     ";
      break;
    case plot:
      break;
    }
    mode = plot;
  }

  enum Mode { command, plot };

  Mode mode;
  std::string filename;
  std::ofstream stream;
  std::string plotDelim;
};
  
#endif // DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH
