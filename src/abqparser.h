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

#ifndef _ABQPARSER_H_
#define _ABQPARSER_H_

#include <stdio.h>
#include <stdlib.h>

#include <vector>
using namespace std;

/*
  A parser for the volumetric mesh text file format.
  Note: the end user never needs to use this class directly.
  See volumetricMesh.h .
*/


class ABQParser
{
public:

    ABQParser();
    ~ABQParser();

    int open(const char * filename);

    // return the next line, s must be externally allocated string
    // if last line, return will be NULL
    char * getNextLine(char * s, int numRetainedSpaces=0);

    void rewindToStart();
    void close();

    static void upperCase(char * s);
    static void stripBlanks(char * s, int numRetainedSpaces); // any whitespace equal in length or longer to "numRetainedSpaces" is shrunk to "numRetainedSpaces" and retained
    static void beautifyLine(char * s, int numRetainedSpaces); // stripBlanks + removes trailing "\n"

protected:

    FILE * fin;
    vector< FILE* > fileStack;
    int currentPos;

    char directoryName[4096];
};



#endif

