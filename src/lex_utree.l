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
%{
#include "parse_utree.h"
#include "mptp.h"

static size_t string_length = 0;

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)realloc(utree_lval.s, (*dstlen + srclen + 1));
  if (!mem)
  {
    fatal("Unable to allocate enough memory.");
  }
  utree_lval.s =  mem;
  strncpy(&utree_lval.s[*dstlen], src, srclen);
  utree_lval.s[*dstlen + srclen] = 0;
  *dstlen += srclen;
  return utree_lval.s;
}

%}
%option noyywrap
%option prefix="utree_"
%option nounput
%option noinput
%s apos
%s quot

%%
<quot>{
\\\"                                         append(&string_length, "\"", 1);
\'                                           append(&string_length, "\'", 1);
\"                                           {BEGIN(INITIAL);return STRING;}
}

<apos>{
\\\'                                         append(&string_length, "\'", 1);
\"                                           append(&string_length, "\"", 1);
\'                                           {BEGIN(INITIAL);return STRING;}
}

<apos,quot>{
\\n                                          append(&string_length, "\n", 1);
\\t                                          append(&string_length, "\t", 1);
\\\\                                         append(&string_length, "\\", 1);
([^\"\'\\]|\n)+                              append(&string_length, utree_text, utree_leng);
}

\:                                           return COLON;
\;                                           return SEMICOLON;
\)                                           return CPAR;
\(                                           return OPAR;
\,                                           return COMMA;
\"                                           {string_length =0; BEGIN(quot);}
\'                                           {string_length =0; BEGIN(apos);}
[\+|\-]?[0-9]+                               {utree_lval.d = xstrndup(utree_text, utree_leng); return NUMBER;}
[\+|\-]?[0-9]+\.[0-9]+([e|E][\+|\-]?[0-9]+)? {utree_lval.d = xstrndup(utree_text, utree_leng); return NUMBER;}
[A-Za-z0-9\'\.\-\_]*                         {utree_lval.s = xstrndup(utree_text, utree_leng); return STRING;}
[ \t\n\r]                                    ;
.                                            {fatal("Syntax error (%c)\n", utree_text[0]);}
%%
