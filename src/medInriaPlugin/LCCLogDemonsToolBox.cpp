/*=========================================================================

 medInria

 Copyright (c) INRIA 2013. All rights reserved.
 See LICENSE.txt for details.
 
  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

=========================================================================*/

#include "LCCLogDemons.h"
#include "LCCLogDemonsToolBox.h"

#include <QtGui>

#include <dtkCore/dtkAbstractDataFactory.h>
#include <dtkCore/dtkAbstractData.h>
#include <dtkCore/dtkAbstractProcessFactory.h>
#include <dtkCore/dtkAbstractProcess.h>
#include <dtkCore/dtkAbstractViewFactory.h>
#include <dtkCore/dtkSmartPointer.h>

#include <medAbstractView.h>
#include <medRunnableProcess.h>
#include <medJobManager.h>

#include <medAbstractDataImage.h>

#include <medToolBoxFactory.h>
#include <medRegistrationSelectorToolBox.h>
#include <medProgressionStack.h>
#include <medGroupBox.h>

class LCCLogDemonsToolBoxPrivate
{
public:
    QComboBox * simMeasureComboBox;
    QSpinBox * stepSizeZ;
    QDoubleSpinBox * minimalVariance;
    QCheckBox * doubleIterations;

    medProgressionStack * progression_stack;
};

LCCLogDemonsToolBox::LCCLogDemonsToolBox(QWidget *parent) : medRegistrationAbstractToolBox(parent), d(new LCCLogDemonsToolBoxPrivate)
{
    this->setTitle("LCCLogDemons");
    
    QWidget *widget = new QWidget(this);
    
    QFormLayout *layout = new QFormLayout(widget);
    
    medGroupBox * advancedBox = new medGroupBox(this);
    advancedBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    advancedBox->setTitle(tr("Advanced Parameters"));
    advancedBox->setCollapsible(true);
    advancedBox->setCheckable(true);
    advancedBox->setChecked(false);
    advancedBox->setStyleSheet("medGroupBox {margin: 10px;}");
    QWidget * advancedWidgets = new QWidget(advancedBox);
    QFormLayout *advancedformLayout = new QFormLayout(advancedWidgets);
    QVBoxLayout * vlayout = new QVBoxLayout(advancedBox);
    vlayout->addWidget(advancedWidgets);

    // Standard parameters
    d->doubleIterations = new QCheckBox();
    layout->addRow(new QLabel(tr("Double "),this),d->doubleIterations);


    // Advanced parameters
    d->minimalVariance = new QDoubleSpinBox();
    advancedformLayout->addRow(new QLabel(tr("Minimal Variance"),this),d->minimalVariance);

    QPushButton *runButton = new QPushButton(tr("Run"), this);

    d->progression_stack = new medProgressionStack(widget);
    QHBoxLayout *progressStackLayout = new QHBoxLayout;
    progressStackLayout->addWidget(d->progression_stack);
    
    this->addWidget(widget);
    this->addWidget(advancedBox);
    this->addWidget(runButton);
    this->addWidget(d->progression_stack);
    
    connect(runButton, SIGNAL(clicked()), this, SLOT(run()));
}

LCCLogDemonsToolBox::~LCCLogDemonsToolBox()
{
    delete d;
    
    d = NULL;
}

bool LCCLogDemonsToolBox::registered()
{
    return medToolBoxFactory::instance()->
    registerToolBox<LCCLogDemonsToolBox>("LCCLogDemonsToolBox",
                               tr("LCC Log Demons"),
                               tr("short tooltip description"),
                               QStringList() << "registration");
}

void LCCLogDemonsToolBox::run()
{
    
    if(!this->parentToolBox())
        return;
    medRegistrationSelectorToolBox * parentTB = this->parentToolBox();
    dtkSmartPointer <dtkAbstractProcess> process;
    
    if (this->parentToolBox()->process() &&
        (parentTB->process()->identifier() == "LCCLogDemons"))
    {
        process = parentTB->process();
        
    }
    else
    {
        process = dtkAbstractProcessFactory::instance()->createSmartPointer("LCCLogDemons");
        parentTB->setProcess(process);
    }
    dtkAbstractData *fixedData = parentTB->fixedData();
    dtkAbstractData *movingData = parentTB->movingData();
    
    
    if (!fixedData || !movingData)
        return;
    
    
    LCCLogDemons *process_Registration = dynamic_cast<LCCLogDemons *>(process.data());
    if (!process_Registration)
    {
        qWarning() << "registration process doesn't exist" ;
        return;
    }
    process_Registration->setUpdateRule(2);
    process_Registration->setVerbosity(true);
    std::vector<unsigned int> iterations;
    iterations.push_back(2); iterations.push_back(0); iterations.push_back(0);
    process_Registration->setIterations(iterations);
    
    process->setInput(fixedData,  0);
    process->setInput(movingData, 1);
    
    medRunnableProcess *runProcess = new medRunnableProcess;
    runProcess->setProcess (process);
    
    d->progression_stack->addJobItem(runProcess, tr("Progress:"));
    //If there is no observer to track the progression,
    //make the progress bar spin:
    //d->progression_stack->setActive(runProcess,true);
    
    connect (runProcess, SIGNAL (success  (QObject*)),  this, SIGNAL (success ()));
    connect (runProcess, SIGNAL (failure  (QObject*)),  this, SIGNAL (failure ()));
    connect (runProcess, SIGNAL (cancelled (QObject*)), this, SIGNAL (failure ()));
    //First have the moving progress bar,
    //and then display the remaining % when known
    connect (runProcess, SIGNAL(activate(QObject*,bool)),
             d->progression_stack, SLOT(setActive(QObject*,bool)));
    
    medJobManager::instance()->registerJobItem(runProcess);
    QThreadPool::globalInstance()->start(dynamic_cast<QRunnable*>(runProcess));
}
