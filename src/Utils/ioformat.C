
#include <sstream>

#include "ioformat.h"

using namespace std;

namespace OPTPP {

ostream& operator << (ostream& ut, oformatstate const& fmt)
{
  ut.width     (fmt.owidth);
  ut.precision (fmt.oprecision);
  ut.fill      (fmt.ofill);
  ut.flags     (fmt.oflags);
  return ut;
}

ostream& operator >> (ostream& ut, oformatstate & fmt)
{
  fmt.owidth     = ut.width();
  fmt.oprecision = ut.precision();
  fmt.ofill      = ut.fill();
  fmt.oflags     = ut.flags();
  return ut;
}

oformatstate::oformatstate (ostream& ut)
  :owidth(ut.width()),
   oprecision(ut.precision()),
   ofill(ut.fill()),
   oflags(ut.flags())
{}

oformatstate::oformatstate (char code, int w, int p, char c, ios::fmtflags f)
  :owidth(w),
   oprecision(p),
   ofill(c),
   oflags(f)
{
  if (owidth < 0)
    {
      oflags |= ios::left;
      owidth = -owidth;
    }
  switch(code)
    {
    case 'd':
    case 'i':
    case 'p':
    case 'u':
    case 'C':
    case 'c':
    case 'S':
    case 's':
      oflags |= ios::dec;
      break;
    case 'o':
      oflags |= ios::oct;
      break;
    case 'X': oflags |= ios::uppercase;
    case 'x':
      oflags |= ios::hex;
      break;
    case 'E': oflags |= ios::uppercase;
    case 'e':
      oflags |= ios::scientific;
      break;
    case 'f':
      oflags |= ios::fixed;
      break;
    case 'G': oflags |= ios::uppercase;
    case 'g':
      break;
	
    }
}

string format(double val, oformatstate const& fmt)
{
  std::ostringstream ut;

  ut << fmt << val;
  std::string s = ut.str(); // copy into a string
  return string(s);
}
string format(int val, oformatstate const& fmt)
{
  std::ostringstream ut;
  
  ut << fmt << val;
  std::string s = ut.str(); // copy into a string
  return string(s);
}

} // namespace OPTPP
