
void npsol_setup(char* string)
{
  int strlen = 72;
#ifndef __xlc__
  npoptn_(string, strlen);
#else
  npoptn(string, strlen);
#endif
}

