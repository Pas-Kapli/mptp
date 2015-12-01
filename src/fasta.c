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
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "mptp.h"

#define MEMCHUNK 4096

/* please note that these functions will return a pointer to a buffer
   allocated here for the query header and sequence. This buffers will
   be overwritten on the next call of query_getnext. */

pll_fasta_t * pll_fasta_open(const char * filename, const unsigned int * map)
{
  int i;
  pll_fasta_t * fd = (pll_fasta_t *)malloc(sizeof(pll_fasta_t));
  if (!fd) return NULL;

  /* allocate space */

  fd->lineno = 0;

  fd->no = -1;

  fd->chrstatus = map;

  /* open file */
  fd->fp = fopen(filename, "r");
  if (!(fd->fp))
  {
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(errmsg, 200, "Unable to open file (%s)", filename);
    free(fd);
    return PLL_FAILURE;
  }

  /* get filesize */
  if (fseek(fd->fp, 0, SEEK_END))
  {
    pll_errno = PLL_ERROR_FILE_SEEK;
    snprintf(errmsg, 200, "Unable to seek in file (%s)", filename);
    free(fd);
    return PLL_FAILURE;
  }
  fd->filesize = ftell(fd->fp);

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(i=0; i<256; i++)
    fd->stripped[i] = 0;

  fd->line[0] = 0;
  if (!fgets(fd->line, PLL_LINEALLOC, fd->fp))
  {
    pll_errno = PLL_ERROR_FILE_SEEK;
    snprintf(errmsg, 200, "Unable to read file (%s)", filename);
    free(fd);
    return PLL_FAILURE;
  }
  fd->lineno = 1;

  return fd;
}

int pll_fasta_rewind(pll_fasta_t * fd)
{
  int i;

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(i=0; i<256; i++)
    fd->stripped[i] = 0;

  fd->line[0] = 0;
  if (!fgets(fd->line, PLL_LINEALLOC, fd->fp))
  {
    pll_errno = PLL_ERROR_FILE_SEEK;
    snprintf(errmsg, 200, "Unable to rewind and cache data");
    return PLL_FAILURE;
  }
  fd->lineno = 1;

  return PLL_SUCCESS;
}

void pll_fasta_close(pll_fasta_t * fd)
{
  fclose(fd->fp);
  free(fd);
}

int pll_fasta_getnext(pll_fasta_t * fd, char ** head,
                                 long * head_len, char ** seq, 
                                 long * seq_len, long * seqno)
{
  void * mem;
  long head_alloc = MEMCHUNK;
  long seq_alloc = MEMCHUNK;
  
  *head_len = 0;
  *seq_len = 0;

  /* allocate sequence buffers */
  *head = (char *)malloc((size_t)(head_alloc));
  if (!(*head))
    return PLL_FAILURE;

  *seq = (char *)malloc((size_t)(seq_alloc));
  if (!(*seq))
  {
    free(*head);
    return PLL_FAILURE;
  }

  /* read line and increase line number */

  while (fd->line[0])
    {
      /* read header */

      if (fd->line[0] != '>')
      {
        pll_errno = PLL_ERROR_FASTA_INVALIDHEADER;
        snprintf(errmsg, 200, "Illegal header line in query fasta file");
        free(*head);
        free(*seq);
        return PLL_FAILURE;
      }

       
      long headerlen;
      if (strchr(fd->line+1,'\r'))
        headerlen = xstrchrnul(fd->line+1, '\r') - (fd->line+1);
      else
        headerlen = xstrchrnul(fd->line+1, '\n') - (fd->line+1);

      *head_len = headerlen;


      if (headerlen + 1 > head_alloc)
      {
        head_alloc = headerlen + 1;
        mem = realloc(*head, (size_t)(head_alloc));
        if (!mem)
        {
          pll_errno = PLL_ERROR_MEM_ALLOC;
          snprintf(errmsg, 200, "Unable to allocate enough memory.");
          free(*head);
          free(*seq);
          return PLL_FAILURE;
        }
        *head = (char *)mem;
      }

      memcpy(*head, fd->line + 1, (size_t)headerlen);
      *(*head + headerlen) = 0;

      /* get next line */

      fd->line[0] = 0;
      if (!fgets(fd->line, PLL_LINEALLOC, fd->fp))
      {
        /* do nothing */
      }
      fd->lineno++;

      /* read sequence */

      *seq_len = 0;

      while (fd->line[0] && (fd->line[0] != '>'))
        {
          char c;
          char m;
          char * p = fd->line;

          while((c = *p++))
            {
              m = (char) fd->chrstatus[(int)c];
              switch(m)
                {
                case 0:
                  /* character to be stripped */
                  fd->stripped_count++;
                  fd->stripped[(int)c]++;
                  break;

                case 1:
                  /* legal character */
                  if (*seq_len + 1 > seq_alloc)
                    {
                      seq_alloc += MEMCHUNK;
                      mem = realloc(*seq, (size_t)(seq_alloc));
                      if (!mem)
                      {
                        pll_errno = PLL_ERROR_MEM_ALLOC;
                        snprintf(errmsg, 200,
                                 "Unable to allocate enough memory.");
                        free(*head);
                        free(*seq);
                        return PLL_FAILURE;
                      }
                      *seq = (char *)mem;
                    }
                  *(*seq + *seq_len) = c;
                  (*seq_len)++;

                  break;

                case 2:
                  /* fatal character */
                  if (c>=32)
                  {
                    pll_errno = PLL_ERROR_FASTA_ILLEGALCHAR;
                    snprintf(errmsg, 200, "illegal character '%c' "
                                              "on line %ld in the fasta file",
                                              c, fd->lineno);
                  }
                  else
                  {
                    pll_errno = PLL_ERROR_FASTA_UNPRINTABLECHAR;
                    snprintf(errmsg, 200, "illegal unprintable character "
                                              "%#.2x (hexadecimal) on line %ld "
                                              "in the fasta file",
                                              c, fd->lineno);
                  }
                  return PLL_FAILURE;

                case 3:
                  /* silently stripped chars */
                  break;

                }
            }

          fd->line[0] = 0;
          if (!fgets(fd->line, PLL_LINEALLOC, fd->fp))
          {
            /* do nothing */
          }
          fd->lineno++;
        }

      /* add zero after sequence */

      if (*seq_len + 1 > seq_alloc)
        {
          seq_alloc += MEMCHUNK;
          mem = realloc(*seq, (size_t)seq_alloc);
          if (!mem)
          {
            pll_errno = PLL_ERROR_MEM_ALLOC;
            snprintf(errmsg, 200, "Unable to allocate enough memory.");
            free(*head);
            free(*seq);
            return PLL_FAILURE;
          }
          *seq = (char *)mem;
        }
      *(*seq + *seq_len) = 0;

      fd->no++;
      *seqno = fd->no;

      return PLL_SUCCESS;
    }

  
  snprintf(errmsg, 200, "End of file\n");
  pll_errno = PLL_ERROR_FILE_EOF;
  free(*head);
  free(*seq);
  return PLL_FAILURE;
}

long pll_fasta_getfilesize(pll_fasta_t * fd)
{
  return fd->filesize;
}

long pll_fasta_getfilepos(pll_fasta_t * fd)
{
  return ftell(fd->fp);
}
