/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "mptp.h"

#define MPTP_RAND48_MULT_0   (0xe66d)
#define MPTP_RAND48_MULT_1   (0xdeec)
#define MPTP_RAND48_MULT_2   (0x0005)
#define MPTP_RAND48_ADD      (0x000b)

static const char * progress_prompt;
static unsigned long progress_next;
static unsigned long progress_size;
static unsigned long progress_chunk;
static const unsigned long progress_granularity = 200;

void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(1);
}

void progress_init(const char * prompt, unsigned long size)
{
  if (!opt_quiet)
  {
    progress_prompt = prompt;
    progress_size = size;
    progress_chunk = size < progress_granularity ?
      1 : size  / progress_granularity;
    progress_next = 0;
    fprintf(stderr, "%s %.0f%%", prompt, 0.0);
  }
}

void progress_update(unsigned int progress)
{
  if (!opt_quiet)
  {
    if (progress >= progress_next)
    {
      fprintf(stderr, "  \r%s %.0f%%", progress_prompt,
              100.0 * progress  / progress_size);
      progress_next = progress + progress_chunk;
    }
  }
}

void progress_done()
{
  if (!opt_quiet)
    fprintf(stderr, "  \r%s %.0f%%\n", progress_prompt, 100.0);
}

#if 0
void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t = NULL;
  if (posix_memalign(& t, alignment, size) == -1)
    fatal("Unable to allocate enough memory.");

  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}
#else
void * xmalloc(size_t size)
{
  void * t;
  t = malloc(size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}
#endif

void * xcalloc(size_t nmemb, size_t size)
{
  void * t;
  t = calloc(nmemb,size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

char * xstrchrnul(char *s, int c)
{
  char * r = strchr(s, c);

  if (r)
    return r;
  else
    return (char *)s + strlen(s);
}

char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
}

char * xstrndup(const char * s, size_t len)
{
  char * p = (char *)xmalloc(len+1);
  strncpy(p,s,len);
  p[len] = 0;
  return p;
}

long getusec(void)
{
  struct timeval tv;
  if(gettimeofday(&tv,0) != 0) return 0;
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

FILE * xopen(const char * filename, const char * mode)
{
  FILE * out = fopen(filename, mode);
  if (!out)
    fatal("Cannot open file %s", filename);

  return out;
}

void random_init(unsigned short * rstate, long seedval)
{
  /* emulate drand48() */
  rstate[0] = 0x330e;
  rstate[1] = seedval & 0xffffl;
  rstate[2] = seedval >> 16;
}

static void mptp_randomize(unsigned short * rstate) 
{
  unsigned long accu;
  unsigned short temp[2];

  accu = MPTP_RAND48_MULT_0 * (unsigned long)rstate[0] +
         MPTP_RAND48_ADD;
  temp[0] = (unsigned short)accu;        /* lower 16 bits */

  accu >>= sizeof(unsigned short) * 8;
  accu += (unsigned long)MPTP_RAND48_MULT_0 * (unsigned long)rstate[1] +
          (unsigned long)MPTP_RAND48_MULT_1 * (unsigned long)rstate[0];
  temp[1] = (unsigned short)accu;        /* middle 16 bits */

  accu >>= sizeof(unsigned short) * 8;
  accu += MPTP_RAND48_MULT_0 * rstate[2] +
          MPTP_RAND48_MULT_1 * rstate[1] +
          MPTP_RAND48_MULT_2 * rstate[0];

  rstate[0] = temp[0];
  rstate[1] = temp[1];
  rstate[2] = (unsigned short)accu;
}

double mptp_erand48(unsigned short * rstate)
{
  mptp_randomize(rstate);

  return ldexp((double)rstate[0], -48) +
         ldexp((double)rstate[1], -32) +
         ldexp((double)rstate[2], -16);
}

long mptp_nrand48(unsigned short * rstate) 
{
  mptp_randomize(rstate);
  return ((long)rstate[2] << 15) + ((long)rstate[1] >> 1);
}
