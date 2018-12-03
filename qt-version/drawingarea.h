#ifndef DRAWINGAREA_H
#define DRAWINGAREA_H

#include <QObject>
#include <QWidget>
#include <QPainter>
#include <QMouseEvent>

#include <cmath>

#include "simulation.h"

class DrawingArea : public QWidget {
    Q_OBJECT

public:
    DrawingArea(QWidget *parent);

    ~DrawingArea();

    void setState(Simulation *simulation);

private:

    Simulation *m_s;

    bool m_changetype;

    QPoint m_selectionBoxStart, m_selectionBoxEnd;
    bool m_selection;

    QPolygon polygon;

protected:
    void paintEvent(QPaintEvent *event);

    void mousePressEvent(QMouseEvent *e);

    void mouseReleaseEvent(QMouseEvent *e);

    void mouseMoveEvent(QMouseEvent *e);

    void mouseDoubleClickEvent(QMouseEvent *e);

    void resizeEvent(QResizeEvent *e);

    void keyPressEvent(QKeyEvent *e);

    signals:
            void

    handChanged();
};

#endif // DRAWINGAREA_H
