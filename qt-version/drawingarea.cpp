#include "drawingarea.h"

DrawingArea::DrawingArea(QWidget *parent) : QWidget(parent) {
    //sthis->setStyleSheet("background-color: white; ");
    m_s = 0;

}

DrawingArea::~DrawingArea() {

}


void DrawingArea::setState(Simulation *simulation) {

    m_s = simulation;
    update();
}

void DrawingArea::paintEvent(QPaintEvent *event) {
    if (m_s == 0 || m_s->display.nw == 0 || m_s->display.nh == 0) return;

    QPainter painter(this);
    //painter.setBackground(QBrush(Qt::green));
    //painter.setBackgroundMode(Qt::OpaqueMode);
    //painter.setBackground(QBrush(Qt::green));
    painter.fillRect(rect(), QBrush(Qt::white));

    painter.setRenderHints(QPainter::SmoothPixmapTransform, true);
    painter.setRenderHints(QPainter::Antialiasing, true);

    int w = width();
    int h = height();
    int n1 = m_s->display.nw;
    int n2 = m_s->display.nh;

    int a1p = w / n1;
    int a2p = h / (std::sqrt(3.0) / 2.0 * n2);

    int a = std::min(a1p, a2p);

    double dx = 0.0;
    double dy = 0.0;
    double lx = a * n1;

    double ly = a * std::sqrt(3.0) / 2.0 * n2;

    dx = (w - lx) / 2.0;
    dy = (h - ly) / 2.0;




    //if (ly>h) a = h/(std::sqrt(3.0)/2.0*n2) ;
    int k = 0;
    for (int i2 = 0; i2 < n2; i2++)
        for (int i1 = 0; i1 < n1; i1++) {

            assert(k < m_s->control.co.rows());

            double x = m_s->control.co(k, 0) * a;
            double y = m_s->control.co(k, 1) * a;

            y = ly - y - a / 2;

            if (m_changetype) {
                int xxmi, yymi, xxma, yyma;
                xxmi = std::min(m_selectionBoxStart.x(), m_selectionBoxEnd.x());
                xxma = std::max(m_selectionBoxStart.x(), m_selectionBoxEnd.x());
                yymi = std::min(m_selectionBoxStart.y(), m_selectionBoxEnd.y());
                yyma = std::max(m_selectionBoxStart.y(), m_selectionBoxEnd.y());

                if (x + dx < xxma && x + dx > xxmi && y + dy < yyma && y + dy > yymi) {
                    m_s->display.changeThisIfShould(k);
                }
            }

            const QColor &color_site = m_s->display.getSiteColor(k);
            const QColor &color_type = m_s->display.getTypeColor(k);

            painter.setPen(QPen(color_site));
            painter.setBrush(QBrush(color_type));
            painter.drawEllipse(x + dx, y + dy, a / 2, a / 2);


            if (m_selection) {
                int xxmi, yymi, xxma, yyma;
                xxmi = std::min(m_selectionBoxStart.x(), m_selectionBoxEnd.x());
                xxma = std::max(m_selectionBoxStart.x(), m_selectionBoxEnd.x());
                yymi = std::min(m_selectionBoxStart.y(), m_selectionBoxEnd.y());
                yyma = std::max(m_selectionBoxStart.y(), m_selectionBoxEnd.y());

                if (x + dx < xxma && x + dx > xxmi && y + dy < yyma && y + dy > yymi) {
                    painter.setPen(QPen(m_s->display.getSiteEditColor()));
                    painter.setBrush(QBrush(m_s->display.getTypeEditColor()));
                    painter.drawEllipse(x + dx - 0.05 * a / 2, y + dy - 0.05 * a / 2, 1.1 * a / 2, 1.1 * a / 2);


                }
            }


            k++;
        }

    if (m_changetype) emit handChanged();
    m_changetype = false;

    painter.setBrush(QBrush());
    painter.setPen(QPen("black"));

    if (m_selection) {
        painter.drawRect(QRectF(m_selectionBoxStart, m_selectionBoxEnd));
    }

    painter.setBrush(QBrush());
    painter.setPen(QPen("green"));
    painter.drawPolygon(polygon);

}

void DrawingArea::mouseDoubleClickEvent(QMouseEvent *e) {
    if (e->button() == Qt::RightButton) {
        if (!polygon.isEmpty()) polygon.pop_back();
        update();
    }

    if (e->button() == Qt::LeftButton) {
        polygon.append(e->pos());
        update();
    }
}

void DrawingArea::resizeEvent(QResizeEvent *e) {
    polygon.clear();
}

void DrawingArea::mousePressEvent(QMouseEvent *e) {
    if (!m_selection) {
        m_selectionBoxStart = e->pos();
        m_selectionBoxEnd = e->pos();
        m_selection = true;
        update();
    }


}

void DrawingArea::mouseReleaseEvent(QMouseEvent *e) {
    m_selectionBoxEnd = e->pos();
    m_selection = false;
    m_changetype = true;
    update();
}

void DrawingArea::mouseMoveEvent(QMouseEvent *e) {
    if (m_selection) {
        m_selectionBoxEnd = e->pos();
        update();
    }
}

void DrawingArea::keyPressEvent(QKeyEvent *e) {

    if (e->key() == Qt::Key_P) {
        if (m_s == 0 || m_s->display.nw == 0 || m_s->display.nh == 0) return;

        int w = width();
        int h = height();
        int n1 = m_s->display.nw;
        int n2 = m_s->display.nh;

        int a1p = w / n1;
        int a2p = h / (std::sqrt(3.0) / 2.0 * n2);

        int a = std::min(a1p, a2p);

        double dx = 0.0;
        double dy = 0.0;
        double lx = a * n1;

        double ly = a * std::sqrt(3.0) / 2.0 * n2;

        dx = (w - lx) / 2.0;
        dy = (h - ly) / 2.0;




        //if (ly>h) a = h/(std::sqrt(3.0)/2.0*n2) ;
        int k = 0;
        for (int i2 = 0; i2 < n2; i2++)
            for (int i1 = 0; i1 < n1; i1++) {

                assert(k < m_s->control.co.rows());

                double x = m_s->control.co(k, 0) * a;
                double y = m_s->control.co(k, 1) * a;
                y = ly - y - a / 2;

                if (polygon.containsPoint(QPoint(x + dx + a / 4, y + dy + a / 4), Qt::OddEvenFill)) {
                    m_s->display.changeThisIfShould(k);
                }
                k++;
            }

        emit handChanged();
        polygon.clear();
        update();
    }
}
