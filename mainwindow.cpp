#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>

#include <iostream>


#include "SystemParams.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    animTimer = new QTimer(this); 
    connect(animTimer, SIGNAL(timeout()), this, SLOT(AnimationThread()));
    connect(ui->runButton,	 SIGNAL(clicked()), this, SLOT(SimulationStart()));
    connect(ui->actionOpen,	 SIGNAL(triggered()), this, SLOT(FileOpen()));

    connect(ui->show_points_cb,	 SIGNAL(stateChanged(int)), this, SLOT(CheckBoxesTriggered()));
    connect(ui->show_image_cb,	 SIGNAL(stateChanged(int)), this, SLOT(CheckBoxesTriggered()));
    connect(ui->show_mag_cb,	 SIGNAL(stateChanged(int)), this, SLOT(CheckBoxesTriggered()));

    ui->show_image_cb->setDisabled(true);
    ui->show_mag_cb->setDisabled(true);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete animTimer;
}

void MainWindow::CheckBoxesTriggered()
{
    SystemParams::show_points = ui->show_points_cb->isChecked();
    SystemParams::show_image = ui->show_image_cb->isChecked();
    SystemParams::show_mag = ui->show_mag_cb->isChecked();
    ui->widget->GetGLWidget()->repaint();
}

void MainWindow::FileOpen()
{
    QString qFilename = QFileDialog::getOpenFileName(this, "/home/azer/Desktop/");
    if(qFilename.isEmpty()) return;

    ui->show_image_cb->setDisabled(false);
    ui->show_mag_cb->setDisabled(false);

    ui->widget->GetGLWidget()->SetImage(qFilename);
    ui->widget->GetGLWidget()->repaint();
    ui->widget->SetScrolls();
}

void MainWindow::AnimationThread()
{
    float duration = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;

    this->ui->widget->GetGLWidget()->SetRunningTime(duration);
    ui->timeLabel->setText("Time: " + QString::number(duration));
    ui->pointsLabel->setText("# Points: " + QString::number(this->ui->widget->GetGLWidget()->GetNPoints()));
    ui->iterLabel->setText("Iteration: " + QString::number(this->ui->widget->GetGLWidget()->GetCurrentIter()));

    this->ui->widget->GetGLWidget()->repaint();
    if(this->ui->widget->GetGLWidget()->IsCalculationDone())
    {
        animTimer->stop();
    }
}

void MainWindow::SimulationStart()
{
    std::cout << "simulation start\n";
    this->ui->widget->GetGLWidget()->StartEvolution();
    this->ui->widget->GetGLWidget()->repaint();

    startTime = std::clock();

    SystemParams::f_b = ui->fb_sb->value();
    SystemParams::f_f = ui->ff_sb->value();
    SystemParams::f_a = ui->fa_sb->value();
    SystemParams::delta_bound = ui->delta_bound_sb->value();

    SystemParams::show_points = ui->show_points_cb->isChecked();
    SystemParams::show_image = ui->show_image_cb->isChecked();

    animTimer->start(ui->delay_sb->value());
}
