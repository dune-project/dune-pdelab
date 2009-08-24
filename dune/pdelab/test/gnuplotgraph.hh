#ifndef DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH
#define DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

class GnuplotGraph
{
public:
  GnuplotGraph(const std::string &prefix)
    : mode(command)
    , filename(prefix+".gnuplot")
    , stream(filename.c_str())
    , datname_(prefix+".dat")
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
    datstream.close();
    std::ostringstream command;
    command << "gnuplot " << filename;
    std::system(command.str().c_str());
  }

  const std::string &datname() const {
    return datname_;
  }

  std::ostream &dat() {
    if(!datstream.is_open())
      datstream.open(datname_.c_str());
    return datstream;
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
  std::string datname_;
  std::ofstream datstream;
  std::string plotDelim;
};
  
#endif // DUNE_PDELAB_TEST_GNUPLOTGRAPH_HH
