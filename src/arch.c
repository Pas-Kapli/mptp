/*
    Copyright (C) 2014-2015 Tomas Flouri, Torbjorn Rognes, Jeff Epler

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

unsigned long arch_get_memused()
{
#ifdef _WIN32

  PROCESS_MEMORY_COUNTERS pmc;
  GetProcessMemoryInfo(GetCurrentProcess(),
                       &pmc,
                       sizeof(PROCESS_MEMORY_COUNTERS));
  return pmc.PeakWorkingSetSize;

#else

  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);

# ifdef __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return r_usage.ru_maxrss;
# else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  return r_usage.ru_maxrss * 1024;
# endif

#endif
}

unsigned long arch_get_memtotal()
{
#ifdef _WIN32

  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&ms);
  return ms.ullTotalPhys;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  long ram = 0;
  size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, NULL, 0) == -1)
    fatal("Cannot determine amount of RAM");
  return ram;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  long phys_pages = sysconf(_SC_PHYS_PAGES);
  long pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1))
    fatal("Cannot determine amount of RAM");
  return pagesize * phys_pages;

#else

  struct sysinfo si;
  if (sysinfo(&si))
    fatal("Cannot determine amount of RAM");
  return si.totalram * si.mem_unit;

#endif
}

long arch_get_cores()
{
#ifdef _WIN32
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwNumberOfProcessors;
#else
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

void arch_get_user_system_time(double * user_time, double * system_time)
{
  *user_time = 0;
  *system_time = 0;
#ifdef _WIN32
  HANDLE hProcess = GetCurrentProcess();
  FILETIME ftCreation, ftExit, ftKernel, ftUser;
  ULARGE_INTEGER ul;
  GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
  ul.u.HighPart = ftUser.dwHighDateTime;
  ul.u.LowPart = ftUser.dwLowDateTime;
  *user_time = ul.QuadPart * 100.0e-9;
  ul.u.HighPart = ftKernel.dwHighDateTime;
  ul.u.LowPart = ftKernel.dwLowDateTime;
  *system_time = ul.QuadPart * 100.0e-9;
#else
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, & r_usage);
  * user_time = r_usage.ru_utime.tv_sec * 1.0 
    + r_usage.ru_utime.tv_usec * 1.0e-6;
  * system_time = r_usage.ru_stime.tv_sec * 1.0 
    + r_usage.ru_stime.tv_usec * 1.0e-6;
#endif
}

void arch_srandom()
{
  /* initialize pseudo-random number generator */
  unsigned int seed = opt_seed;
  if (seed == 0)
    {
#ifdef _WIN32
      srand(GetTickCount());
#else
      int fd = open("/dev/urandom", O_RDONLY);
      if (fd < 0)
        fatal("Unable to open /dev/urandom");
      if (read(fd, & seed, sizeof(seed)) < 0)
        fatal("Unable to read from /dev/urandom");
      close(fd);
      srandom(seed);
#endif
    }
  else
    {
#ifdef _WIN32
      srand(seed);
#else
      srandom(seed);
#endif
    }
}

long arch_random()
{
#ifdef _WIN32
  return rand();
#else
  return random();
#endif
}
