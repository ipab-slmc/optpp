// Uncomment this only is you are REALLY sure templates work on your system
//#define WANT_TEMPLATES  

#ifndef ioformat_h
#define ioformat_h

//Convert to string 10/05/01

#include <iostream>
#include <string>

using std::ostream;
using std::string;

#if (defined(__GNUC__) && __GNUC__ >= 3)
typedef std::ios_base::fmtflags   opt_mode;
#else
typedef std::ios::fmtflags opt_mode;
#endif


namespace OPTPP {

class oformatstate
{
public:
  oformatstate(ostream& ut);
#if ( defined(__GNUC__) && __GNUC__ >= 3 || defined (__ICC) )
  oformatstate (char code, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed);
#else
  oformatstate (char code, int w=0, int p=0, char c=' ', opt_mode f=0);
#endif
  friend ostream& operator << (ostream& ut, oformatstate const& fmt);
  friend ostream& operator >> (ostream& ut, oformatstate& fmt);
  
private:
  int owidth;
  int oprecision;
  char ofill;
  opt_mode oflags;
};

string format(double val, oformatstate const& fmt);
string format(int    val, oformatstate const& fmt);

inline string
#if ( defined(__GNUC__) && __GNUC__ >= 3 || defined (__ICC) )
d(int val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
#else
d(int val, int w=0, int p=0, char c=' ', opt_mode f=0)
#endif
{
  oformatstate fmt('d', w, p, c, f);
  return format(val, fmt);
}
inline string
#if ( defined(__GNUC__) && __GNUC__ >= 3 || defined (__ICC) )
e(double val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
#else
e(double val, int w=0, int p=0, char c=' ', opt_mode f=0)
#endif
{
  oformatstate fmt('e', w, p, c, f);
  return format(val, fmt);
}
inline string
#if ( defined(__GNUC__) && __GNUC__ >= 3 || defined (__ICC) )
f(double val, int w=0, int p=0, char c=' ', opt_mode f = std::ios_base::fixed)
#else
f(double val, int w=0, int p=0, char c=' ', opt_mode f=0)
#endif
{
  oformatstate fmt('f', w, p, c, f);
  return format(val, fmt);
}

} // namespace OPTPP

#endif
