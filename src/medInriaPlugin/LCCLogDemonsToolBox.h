/*=========================================================================

 medInria

 Copyright (c) INRIA 2013. All rights reserved.
 See LICENSE.txt for details.
 
  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/

#pragma once

#include <medRegistrationAbstractToolBox.h>
#include "LCCLogDemonsPluginExport.h"

class LCCLogDemonsToolBoxPrivate;

class LCCLogDemonsPLUGIN_EXPORT LCCLogDemonsToolBox : public medRegistrationAbstractToolBox
{
    Q_OBJECT
    
public:
    LCCLogDemonsToolBox(QWidget *parent = 0);
    ~LCCLogDemonsToolBox();
    
public:
    static bool registered();
    
public slots:
    void chooseUpdateRule(int);
    void run();
    
private:
    LCCLogDemonsToolBoxPrivate *d;
};


