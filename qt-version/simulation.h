
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <list>
#include <set>
#include <cstdarg>
#include <tuple>
#include <stdexcept>
#include <random>

#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <Eigen/Sparse>

#include "growth_base.h"


#ifndef SIMULATION_H
#define SIMULATION_H

struct Display {
    int nw, nh;

    Eigen::VectorXi types;
    Eigen::VectorXi sites;
    bool hasHandChanged;

    QVector <QColor> cm_state;
    QVector <QColor> cm_sites;

    int edit_type;
    int edit_site;

    bool doEdit_type;
    bool doEdit_site;

    QColor edit_type_color;
    QColor edit_site_color;

    inline void updateFrom(const Control &control) {
        nw = control.nw;
        nh = control.nh;
        types = control.types;
        sites = control.sites;
        hasHandChanged = false;
    }

    inline bool hasChangedByHand() const { return hasHandChanged; }

    inline void getColorMap(const Control &control) {
        cm_state.clear();
        for (auto c : control.color_state) cm_state.push_back(c.c_str());
        cm_sites.clear();
        for (auto c : control.color_sites) cm_sites.push_back(c.c_str());
    }

    inline void setEditColors(const std::string &tlab, bool doEditType,
                              const std::string &slab, bool doEditSite,
                              const Parameters &parameters) {
        try {
            edit_type = parameters.lab2type.at(tlab);
            edit_site = parameters.lab2site.at(slab);
            edit_type_color = cm_state[edit_type];
            edit_site_color = cm_sites[edit_site];
        } catch (const std::exception &e) {
            edit_type = 0;
            edit_site = 0;
            edit_type_color = QColor("white");
            edit_site_color = QColor("white");

        }
        doEdit_site = doEditSite;
        doEdit_type = doEditType;
    }

    inline const QColor &getTypeColor(int k) {
        assert(k < types.rows());
        assert(types[k] < cm_state.size());

        return cm_state[types[k]];
    }

    inline const QColor &getSiteColor(int k) {

        assert(k < sites.rows());
        assert(sites[k] < cm_sites.size());

        return cm_sites[sites[k]];
    }

    inline const QColor &getTypeEditColor() {
        return edit_type_color;
    }

    inline const QColor &getSiteEditColor() {
        return edit_site_color;
    }

    inline void changeThisIfShould(int k) {
        assert(k < types.rows());
        assert(k < sites.rows());

        if (doEdit_type) types(k) = edit_type;
        if (doEdit_site) sites(k) = edit_site;
    }
};

class Simulation {
public:
    std::ofstream ofs;

    Label state, site;
    Parameters parameters;
    Control control;

    Display display;


    inline void close() {
        ofs.close();

    }

    inline void open(std::string fn) {
        if (!ofs.is_open()) ofs.open(fn.c_str());
        ofs.flush();
    }

    inline void write() {
        if (ofs.is_open()) control.logStream(ofs, parameters);
    }


public:
    //Simulation();
    //~Simulation()
};

#endif // SIMULATION_H
