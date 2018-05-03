
/*
 * C file utilities needed by PDS
 *
 * WARNING:  These are PDS-specific!!!
 *
 */

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>

#include "pds.h"

int bin_open(char *filename, int *fd)
{
  *fd = open(filename, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  if (*fd == -1)
    return -1;
  else
    return 0;
}

int bin_close(int fd)
{
  int error;

  error = close(fd);
  return error;
}
#ifdef needed
int read_fopen(char *filename, int *fd)
{
  if ((*fd = fopen(filename, "r")) == NULL)
    return 9;
  else
    return 1;
}

int file_fclose(int fd)
{
  fclose(fd);
  return 1;
}

int read_next(char *filename, int pos, int *record, int size, int length, int *end)
{
  fseek(filename, (long) pos, SEEK_SET);
  fread(READ_TYPE record, size, length, filename);
  end = feof(filename);
  return 1;
}
#endif
