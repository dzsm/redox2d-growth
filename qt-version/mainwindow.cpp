
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

void MainWindow::updateSimulationControl() {
    std::ostringstream controlStream;
    controlStream << "time " << ui->lineEdit_time->text().toStdString() << std::endl;
    controlStream << "seed " << ui->lineEdit_seed->text().toStdString() << std::endl;
    controlStream << "run " << ui->lineEdit_run->text().toStdString() << std::endl;
    controlStream << "every " << ui->lineEdit_every->text().toStdString() << std::endl;
    controlStream << "flagPoisson " << ui->checkBox_poisson->isChecked() << std::endl;
    controlStream << "flagDiffusion " << ui->checkBox_diffusion->isChecked() << std::endl;
    controlStream << "flagRedox " << ui->checkBox_redox->isChecked() << std::endl;
    controlStream << "flagPattern " << ui->checkBox_pattern->isChecked() << std::endl;
    controlStream << "colorState " << ui->lineEdit_cmstate->text().toStdString() << std::endl;
    controlStream << "colorSites " << ui->lineEdit_cmsite->text().toStdString() << std::endl;
    controlStream << "poisson_maxiter " << ui->lineEdit_maxiter->text().toStdString() << std::endl;
    controlStream << "poisson_tolerance " << ui->lineEdit_tolerance->text().toStdString() << std::endl;
    controlStream << "poisson_doBICG " << ui->checkBox_mgbicg->isChecked() << std::endl;

    simulation.control.fromString(controlStream.str());

    ui->plainTextEdit_control->clear();
    ui->plainTextEdit_control->appendPlainText(controlStream.str().c_str());

    simulation.display.updateFrom(simulation.control);
    simulation.display.getColorMap(simulation.control);
    simulation.display.setEditColors(ui->lineEdit_selstate->text().trimmed().toStdString(),
                                     ui->checkBox_state->isChecked(),
                                     ui->lineEdit_selsite->text().trimmed().toStdString(),
                                     ui->checkBox_site->isChecked(),
                                     simulation.parameters);
}

void MainWindow::displaySimulationVariables() {
    ui->lineEdit_time->setText(QString::number(simulation.control.variables["t"]));
    ui->lineEdit_var0->setText(QString::number(simulation.control.variables["$0"]));
    ui->lineEdit_var1->setText(QString::number(simulation.control.variables["$1"]));
    ui->lineEdit_var2->setText(QString::number(simulation.control.variables["$2"]));
    ui->lineEdit_var3->setText(QString::number(simulation.control.variables["$3"]));
    ui->lineEdit_var4->setText(QString::number(simulation.control.variables["$4"]));
    ui->lineEdit_var5->setText(QString::number(simulation.control.variables["$5"]));
    ui->lineEdit_var6->setText(QString::number(simulation.control.variables["$6"]));
    ui->lineEdit_var7->setText(QString::number(simulation.control.variables["$7"]));
    ui->lineEdit_var8->setText(QString::number(simulation.control.variables["$8"]));
    ui->lineEdit_var9->setText(QString::number(simulation.control.variables["$9"]));

}

MainWindow::MainWindow(QWidget *parent) :
        QMainWindow(parent),
        ui(new Ui::MainWindow) {
    ui->setupUi(this);

    //ui->lineEdit_maxiter->setText(QString::number(simulation.grid.solver.solver.iterations()));
    //ui->lineEdit_tolerance->setText(QString::number(simulation.grid.solver.solver.tolerance()));

    ////////////
    updateSimulationControl();
    displaySimulationVariables();
    /////////////

    connect(ui->lineEdit_cmsite, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_cmstate, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_every, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_seed, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_run, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_time, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_maxiter, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_tolerance, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_selsite, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->lineEdit_selstate, SIGNAL(editingFinished()), this, SLOT(updateSimulationControl()));
    connect(ui->widget, SIGNAL(handChanged()), this, SLOT(stateHandChanged()));


}

MainWindow::~MainWindow() {
    delete ui;
}


void MainWindow::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_0) {
        for (int loop = 0; loop < 1; loop++) {
            simulation.control.step(simulation.parameters);
        }
        simulation.display.updateFrom(simulation.control);
        ui->widget->setState(&simulation);
        displaySimulationVariables();
    }
    if (event->key() == Qt::Key_1) {
        for (int loop = 0; loop < 10; loop++) {
            simulation.control.step(simulation.parameters);
        }
        simulation.display.updateFrom(simulation.control);
        ui->widget->setState(&simulation);
        displaySimulationVariables();
    }
    if (event->key() == Qt::Key_2) {
        for (int loop = 0; loop < 100; loop++) {
            simulation.control.step(simulation.parameters);
        }
        simulation.display.updateFrom(simulation.control);
        ui->widget->setState(&simulation);
        displaySimulationVariables();
    }
    if (event->key() == Qt::Key_3) {
        for (int loop = 0; loop < 1000; loop++) {
            simulation.control.step(simulation.parameters);
        }
        simulation.display.updateFrom(simulation.control);
        ui->widget->setState(&simulation);
        displaySimulationVariables();
    }
}

void MainWindow::on_pushButton_load_st_clicked() {
    QFile file(ui->lineEdit_fnst->text().toStdString().c_str());
    if (!file.open(QFile::ReadOnly | QFile::Text));

    QTextStream ts(&file);

    ui->plainTextEdit_state->clear();
    ui->plainTextEdit_state->appendPlainText(ts.readAll());
}

void MainWindow::on_pushButton_save_st_clicked() {
    QFile file(ui->lineEdit_fnst->text().toStdString().c_str());
    if (!file.open(QFile::WriteOnly | QFile::Text));

    file.write(ui->plainTextEdit_state->toPlainText().toStdString().c_str());


}

void MainWindow::on_pushButton_get_st_clicked() {
    std::ostringstream os;
    simulation.state.toStream(os, "state");
    simulation.site.toStream(os, "sites");
    ui->plainTextEdit_state->appendPlainText(os.str().c_str());
}

void MainWindow::on_pushButton_set_st_clicked() {
    simulation.state.fromString(ui->plainTextEdit_state->toPlainText().toStdString(), "state");
    simulation.site.fromString(ui->plainTextEdit_state->toPlainText().toStdString(), "sites");
    simulation.control.setup(simulation.state, simulation.site, simulation.parameters);
    simulation.display.updateFrom(simulation.control);
    ui->widget->setState(&simulation);

    displaySimulationVariables();
}


void MainWindow::on_pushButton_load_prm_clicked() {
    QFile file(ui->lineEdit_fnprm->text().toStdString().c_str());
    if (!file.open(QFile::ReadOnly | QFile::Text));

    QTextStream ts(&file);

    ui->plainTextEdit_prm->clear();
    ui->plainTextEdit_prm->appendPlainText(ts.readAll());

}

void MainWindow::on_pushButton_save_prm_clicked() {
    QFile file(ui->lineEdit_fnprm->text().toStdString().c_str());
    if (!file.open(QFile::WriteOnly | QFile::Text));

    QTextStream ts(&file);

    file.write(ui->plainTextEdit_prm->toPlainText().toStdString().c_str());
}


void MainWindow::on_pushButton_set_prm_clicked() {
    simulation.parameters.fromString(ui->plainTextEdit_prm->toPlainText().toStdString());
    simulation.control.setup(simulation.state, simulation.site, simulation.parameters);
    simulation.display.updateFrom(simulation.control);
    ui->widget->setState(&simulation);

    //simulation.parameters.logStream(std::cout);
    displaySimulationVariables();
}



//void MainWindow::on_lineEdit_sel_editingFinished()
//{
//int rt = simulation.grid.site2type[ui->lineEdit_sel->text().toStdString()];
//ui->widget->setState(simulation.grid.n1,simulation.grid.n2,simulation.grid.types,m_cm,rt);

//}

void MainWindow::on_checkBox_save_toggled(bool checked) {
    if (checked) ui->checkBox_open->toggle();

}

void MainWindow::on_checkBox_open_toggled(bool checked) {
    if (checked) {
        simulation.open(ui->lineEdit_outfn->text().toStdString());
    } else {
        simulation.close();
    }
}

void MainWindow::on_pushButton_clicked() {
    for (int loop = 0; loop < simulation.control.run; loop++) {
        simulation.control.step(simulation.parameters);
        if (loop % simulation.control.every == 0) {
            simulation.write();
            simulation.display.updateFrom(simulation.control);
            ui->widget->setState(&simulation);
            displaySimulationVariables();
        }
    }
    simulation.display.updateFrom(simulation.control);
    ui->widget->setState(&simulation);
    displaySimulationVariables();

}

void MainWindow::on_checkBox_poisson_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_checkBox_redox_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_checkBox_diffusion_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_checkBox_pattern_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_checkBox_state_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_checkBox_site_toggled(bool checked) {
    updateSimulationControl();
}

void MainWindow::on_pushButton_get_prm_clicked() {
    std::string str;
    simulation.parameters.toString(str);
    ui->plainTextEdit_prm->clear();
    ui->plainTextEdit_prm->appendPlainText(str.c_str());
}

void MainWindow::stateHandChanged() {
    simulation.state.modifyById(simulation.display.types, simulation.parameters.type2lab);
    simulation.site.modifyById(simulation.display.sites, simulation.parameters.site2lab);
    simulation.control.setup(simulation.state, simulation.site, simulation.parameters);
    simulation.display.updateFrom(simulation.control);
}

void MainWindow::on_pushButton_plotU_clicked() {

    auto customPlot = ui->widget_plotU;


    // configure axis rect:
    customPlot->setInteractions(
            QCP::iRangeDrag | QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");

    // set up the QCPColorMap:
    QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
    customPlot->addPlottable(colorMap);
    int nx = simulation.control.nw;
    int ny = simulation.control.nh;
    colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
    colorMap->data()->setRange(QCPRange(0, 7), QCPRange(0,
                                                        7)); // and span the coordinate range -4..4 in both key (x) and value (y) dimensions
    // now we assign some data, by accessing the QCPColorMapData instance of the color map:
    double x, y, z;
    int k = 0;
    for (int xIndex = 0; xIndex < simulation.control.nw; ++xIndex) {
        for (int yIndex = 0; yIndex < simulation.control.nh; ++yIndex) {
            x = simulation.control.co(k, 0);
            y = simulation.control.co(k, 1);
            int xI, yI;
            colorMap->data()->coordToCell(x, y, &xI, &yI);

            z = std::real(simulation.control.u(k));
            //colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
            //double r = 3*qSqrt(x*x+y*y)+1e-2;
            //z = 2*x*(qCos(r+2)/r-qSin(r+2)/r); // the B field strength of dipole radiation (modulo physical constants)
            colorMap->data()->setCell(xI, yI, z);
            k++;
        }
    }
    colorMap->data()->recalculateDataBounds();

    // add a color scale:
    QCPColorScale *colorScale = new QCPColorScale(customPlot);
    customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
    colorScale->setType(
            QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
    colorMap->setColorScale(colorScale); // associate the color map with the color scale
    colorScale->axis()->setLabel("Magnetic Field Strength");

    // set the color gradient of the color map to one of the presets:
    colorMap->setGradient(QCPColorGradient::gpPolar);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    // rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
    colorMap->rescaleDataRange();

    // make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
    QCPMarginGroup *marginGroup = new QCPMarginGroup(customPlot);
    customPlot->axisRect()->setMarginGroup(QCP::msBottom | QCP::msTop, marginGroup);
    colorScale->setMarginGroup(QCP::msBottom | QCP::msTop, marginGroup);

    // rescale the key (x) and value (y) axes so the whole color map is visible:
    customPlot->rescaleAxes();
}

void MainWindow::on_checkBox_mgbicg_toggled(bool checked) {
    updateSimulationControl();
}
