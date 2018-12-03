#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "simulation.h"
#include <QKeyEvent>
#include <QTextStream>
#include <QFile>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);

    ~MainWindow();

private
    slots:
            void

    on_pushButton_load_st_clicked();

    void on_pushButton_load_prm_clicked();

    void on_pushButton_save_prm_clicked();

    void on_pushButton_save_st_clicked();

    void on_pushButton_set_prm_clicked();

    void on_pushButton_set_st_clicked();

    void on_pushButton_get_st_clicked();

    void on_checkBox_save_toggled(bool checked);

    void on_checkBox_open_toggled(bool checked);

    void on_pushButton_clicked();

    void updateSimulationControl();

    void displaySimulationVariables();

    void on_checkBox_poisson_toggled(bool checked);

    void on_checkBox_redox_toggled(bool checked);

    void on_checkBox_diffusion_toggled(bool checked);

    void on_checkBox_pattern_toggled(bool checked);

    void on_checkBox_state_toggled(bool checked);

    void on_checkBox_site_toggled(bool checked);

    void on_pushButton_get_prm_clicked();

    void stateHandChanged();

    void on_pushButton_plotU_clicked();

    void on_checkBox_mgbicg_toggled(bool checked);

private:
    Ui::MainWindow *ui;

    Simulation simulation;
    QVector <QColor> m_cm;


    void keyPressEvent(QKeyEvent *event);

};

#endif // MAINWINDOW_H
