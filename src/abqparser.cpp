/*************************************************************************
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Doug L. James, Jovan Popovic                 *
 * Funding: NSF, Link Foundation, Singapore-MIT GAMBIT Game Lab          *
 * Version 2.0                                                           *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "abqparser.h"

#include <string.h>

// convert string to uppercase
void ABQParser::upperCase(char * s)
{  
  char diff = 'A' - 'a';
  for(unsigned int i=0; i<strlen(s); i++)
  {
    if ((s[i] >= 'a') && (s[i] <= 'z'))
      s[i] += diff;
  }
}

// removes all whitespace characters from string s
void ABQParser::stripBlanks(char * s, int numRetainedSpaces)
{
  char * w = s;
  while (*w != '\0')
  {
    while (1) 
    {
      bool cond = (*w == ' ');
      for(int i=1; i<=numRetainedSpaces; i++)
      {
        if (*(w+i) == 0)
        {
          cond = false;
          break;
        }
        cond = (cond && ((w == s) || (*(w+i) == ' '))); // always erase spaces at the beginning of line
      }

      if (!cond)
        break;

      // erase blank

      char * u = w;
      while (*u != '\0') // shift everything left one char
      {
        *u = *(u+1);
        u++;
      }
    }
    w++;
  }
}

ABQParser::ABQParser()
{
  fin = NULL;
  currentPos=-1;
  fileStack.empty();
}

int ABQParser::open(const char * filename)
{

  // extract directory name and filename
  // seek for last '/' in filename
  // if no '/', then directory name is "."
  // else, everything before '/' is directory name, and everything after is filename

  char lfilename[1024];
  strcpy(lfilename,filename);
  char * iter = lfilename;

  char * lastPos = NULL;

  while (*iter != '\0')
  {
    if ((*iter == '/') || (*iter == '\\'))
    {
      lastPos = iter;
    }
    iter++;
  }

  // lastPos is NULL if no slash found, otherwise it points to the last slash
  if (lastPos == NULL)
  {
    directoryName[0] = '.';
    directoryName[1] = '\0';
  }
  else
  {
    iter = lfilename;
    int charPos = 0;
    while (iter != lastPos)
    {
      directoryName[charPos] = *iter;
      charPos++;
      iter++;
    }
    directoryName[charPos] = '\0';
  }

  // now, directoryName has been extracted

  fin = fopen(filename,"r");
  if (!fin)
    return 1;

  currentPos=0;
  fileStack.push_back(fin);

  return 0;
}

ABQParser::~ABQParser()
{
  close();
}

void ABQParser::rewindToStart()
{
  if (currentPos < 0) // no files currently opened
	return;

  while (currentPos > 0)
  {
    fclose(fileStack[currentPos]);
    fileStack.pop_back();
    currentPos--;
  }

  // now, we have currentPos == 0
  fin = fileStack[0];
  rewind(fin);
}

void ABQParser::close()
{
  while (currentPos >= 0)
  {
    fclose(fileStack[currentPos]);
    fileStack.pop_back();
    currentPos--;
  }
}

void ABQParser::beautifyLine(char * s, int numRetainedSpaces)
{  
  stripBlanks(s, numRetainedSpaces);
  //upperCase(s);

  // remove trailing '\n'
  char * pos = s;
  while( *pos != '\0' )
    pos++;

  if ( pos != s ) // string is not zero length
  {
    if( *(pos-1) == '\n' )
      *(pos-1) = '\0';
  }
}

char * ABQParser::getNextLine(char * s, int numRetainedSpaces)
{
  char * code;
  while (((code = fgets(s,4096,fin)) == NULL) && (currentPos > 0)) // if EOF, pop previous file from stack
  {
    currentPos--;
    fclose(fin);
    fileStack.pop_back();
    fin = fileStack[currentPos];
  }

  if (code == NULL) // reached end of main file
  {
    return NULL; 
  }

  beautifyLine(s, numRetainedSpaces);

  // handle the case where input is specified via *INPUT
  while(strncmp(s,"*INCLUDE,INPUT=",15)==0)
  {
    // open up the new file
    char * newFile = &(s[15]);
    char newFileCompleteName[4096];
    sprintf(newFileCompleteName,"%s/%s",directoryName,newFile);

    FILE * finNew = fopen(newFileCompleteName,"r");

    if (!finNew)
    {
      printf("Error: couldn't open include file %s.\n",newFileCompleteName);
      exit(1);
    }

    if ((code = fgets(s,4096,finNew)) != NULL) // new file is not empty
    {
      beautifyLine(s, numRetainedSpaces);

      // register the new file
      currentPos++;
      fileStack.push_back(finNew);
      fin = finNew;
    }
    else
    {
      fclose(finNew);
      printf("Warning: include file is empty.\n");
      code = getNextLine(s, numRetainedSpaces);
      beautifyLine(s, numRetainedSpaces);
      return code;
    }
  }

  return code;
}

